"""Train and run XRI's single-class Mask R-CNN for detector spots."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage
import torch
from torch import nn
from torchvision.models.detection import MaskRCNN
from torchvision.models.detection.anchor_utils import AnchorGenerator
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor, TwoMLPHead
from torchvision.models.detection.mask_rcnn import MaskRCNNHeads, MaskRCNNPredictor
from torchvision.ops import MultiScaleRoIAlign


TILE = 192
MODEL_SCHEMA = "xri.roi.maskrcnn.v1"
DETECTORS = ("D_LAMBDA", "D_MINIPIX")


def build_model():
    # A small stride-four backbone keeps small spots visible and allows CPU training.
    layers = []
    for incoming, outgoing, stride in [(3, 16, 2), (16, 24, 2), (24, 32, 1), (32, 32, 1)]:
        layers.extend((nn.Conv2d(incoming, outgoing, 3, stride, 1), nn.GroupNorm(4, outgoing), nn.SiLU()))
    backbone = nn.Sequential(*layers)
    backbone.out_channels = 32
    return MaskRCNN(
        backbone, num_classes=None, min_size=TILE, max_size=TILE,
        image_mean=[0.0] * 3, image_std=[1.0] * 3,
        rpn_anchor_generator=AnchorGenerator(((6, 12, 24, 40, 64, 96),), ((.2, .5, 1., 2., 5.),)),
        rpn_pre_nms_top_n_train=700, rpn_post_nms_top_n_train=220,
        rpn_pre_nms_top_n_test=500, rpn_post_nms_top_n_test=160,
        rpn_batch_size_per_image=128, box_batch_size_per_image=96,
        box_roi_pool=MultiScaleRoIAlign(["0"], 5, 2),
        box_head=TwoMLPHead(32 * 25, 96), box_predictor=FastRCNNPredictor(96, 2),
        mask_roi_pool=MultiScaleRoIAlign(["0"], 10, 2),
        mask_head=MaskRCNNHeads(32, [32, 32, 32], 1),
        mask_predictor=MaskRCNNPredictor(32, 24, 2),
        box_score_thresh=.03, box_detections_per_img=70,
    )


def read_detector(path, detector, frame=-1, start=0, stop=None):
    with h5py.File(path, "r") as source:
        data = source[f"/entry/data/{detector}"]
        if data.ndim != 3 or not all(1 <= size <= 2048 for size in data.shape[1:]):
            raise ValueError("Expected a frame × row × column detector stack of at most 2048 × 2048.")
        if frame >= 0:
            return np.asarray(data[frame], dtype=np.float64)
        end = data.shape[0] if stop is None else min(data.shape[0], stop)
        if not 0 <= start < end:
            raise ValueError("The selected range has no detector frames.")
        image = np.zeros(data.shape[1:], dtype=np.float64)
        for index in range(start, end, 8):
            image += data[index:min(index + 8, end)].sum(axis=0, dtype=np.float64)
        return image


def features(raw):
    values = np.maximum(np.where(np.isfinite(raw), raw, 0), 0).astype(np.float32)
    root = np.sqrt(values)
    local = ndimage.gaussian_filter(root, 8)
    signal = np.maximum(root - .85 * local, 0)
    scale = max(float(np.percentile(signal, 99.8)), 1.)
    global_signal = np.maximum(root - np.median(root), 0)
    return np.stack((np.clip(signal / scale, 0, 1),
                     np.clip(global_signal / max(np.percentile(global_signal, 99.8), 1.), 0, 1),
                     np.clip(np.log1p(values) / max(np.percentile(np.log1p(values), 99.8), 1.), 0, 1)))


def synthetic(rng):
    yy, xx = np.mgrid[:TILE, :TILE]
    background = rng.uniform(.2, 25) + rng.uniform(0, 12) * xx / TILE + rng.uniform(0, 8) * yy / TILE
    background += rng.uniform(0, 10) * np.exp(-((xx - rng.uniform(0, TILE)) / rng.uniform(14, 65))**2)
    masks, centers = [], []
    for _ in range(int(rng.integers(1, 23)) if rng.random() > .05 else 0):
        x, y = rng.uniform(5, TILE - 5, size=2)
        if any((x - a)**2 + (y - b)**2 < 11**2 for a, b in centers):
            continue
        centers.append((x, y))
        minor = rng.uniform(1.3, 5.)
        major = minor if rng.random() < .25 else minor * rng.uniform(1.15, 4.5)
        angle = rng.uniform(0, np.pi)
        c, s = np.cos(angle), np.sin(angle)
        radius = ((c * (xx - x) + s * (yy - y)) / major)**2 + ((-s * (xx - x) + c * (yy - y)) / minor)**2
        amplitude = np.exp(rng.uniform(np.log(18), np.log(2600)))
        background += amplitude * np.exp(-radius / 2)
        masks.append(radius <= rng.uniform(3.2, 6.2))
    image = rng.poisson(background).astype(np.float32)
    for _ in range(int(rng.integers(0, 35))):
        image[rng.integers(TILE), rng.integers(TILE)] += rng.uniform(3, 900)
    return image, np.asarray(masks, dtype=bool).reshape(-1, TILE, TILE)


def training_sample(image, masks):
    boxes, instances = [], []
    for mask in masks:
        y, x = np.nonzero(mask)
        if x.size >= 8 and np.ptp(x) >= 2 and np.ptp(y) >= 2:
            boxes.append([x.min(), y.min(), x.max() + 1, y.max() + 1])
            instances.append(mask)
    return torch.from_numpy(features(image)), {
        "boxes": torch.as_tensor(np.asarray(boxes).reshape(-1, 4), dtype=torch.float32),
        "labels": torch.ones(len(boxes), dtype=torch.int64),
        "masks": torch.as_tensor(np.asarray(instances).reshape(-1, *image.shape), dtype=torch.uint8),
    }


def weak_masks(image, rectangles):
    masks = []
    for x1, y1, x2, y2 in rectangles:
        x1, x2 = np.clip([x1, x2], 0, image.shape[1])
        y1, y2 = np.clip([y1, y2], 0, image.shape[0])
        if x2 - x1 < 3 or y2 - y1 < 3:
            continue
        patch = image[y1:y2, x1:x2]
        background = np.percentile(patch, 30)
        smooth = ndimage.gaussian_filter(patch.astype(np.float64), .85)
        peak = smooth.max()
        if peak - background < max(4 * np.sqrt(max(background, 1.)), 5.):
            continue
        candidate = smooth > background + max(3 * np.sqrt(max(background, 1.)), .14 * (peak - background))
        labels, _ = ndimage.label(candidate)
        selected = labels.flat[smooth.argmax()]
        if not selected:
            continue
        mask = np.zeros(image.shape, dtype=bool)
        mask[y1:y2, x1:x2] = ndimage.binary_fill_holes(labels == selected)
        if 8 <= mask.sum() <= .85 * patch.size:
            masks.append(mask)
    return np.asarray(masks, dtype=bool).reshape(-1, *image.shape)


def real_crops(args):
    if args.real_masks:
        for path in args.real_masks:
            with np.load(path, allow_pickle=False) as data:
                image, masks = data["image"], data["masks"].astype(bool)
            if image.ndim != 2 or masks.ndim != 3 or masks.shape[1:] != image.shape or not masks.any():
                raise ValueError(f"Invalid real mask annotation: {path}")
            yield image, masks
    if args.real_scan:
        if not args.roi_boxes:
            raise ValueError("--real-scan requires --roi-boxes.")
        records = json.loads(Path(args.roi_boxes).read_text(encoding="utf-8-sig"))["crystals"]
        for detector in DETECTORS:
            rectangles = [[r[k] for k in ("x1", "y1", "x2", "y2")] for r in records
                          if r["detector"] == detector and r["state"] == "valid"]
            if not rectangles:
                continue
            with h5py.File(args.real_scan, "r") as source:
                if f"/entry/data/{detector}" not in source:
                    continue
                frames = source[f"/entry/data/{detector}"].shape[0]
            # Sample short blocks spread over the scan and one summed image.
            block = max(1, min(32, frames // 10))
            spans = [(int(i), int(i) + block) for i in np.linspace(0, frames - block, 6)]
            for begin, end in [(0, frames), *spans]:
                image = read_detector(args.real_scan, detector, start=begin, stop=end)
                masks = weak_masks(image, rectangles)
                print(f"{detector} frames {begin}:{end}: {len(masks)} weak mask labels", flush=True)
                if len(masks):
                    yield image, masks


def crop_real(rng, examples):
    image, masks = examples[int(rng.integers(len(examples)))]
    pad_y, pad_x = (max(0, TILE - size) for size in image.shape)
    if pad_y or pad_x:
        image = np.pad(image, ((0, pad_y), (0, pad_x)))
        masks = np.pad(masks, ((0, 0), (0, pad_y), (0, pad_x)))
    mask = masks[int(rng.integers(len(masks)))]
    ys, xs = np.nonzero(mask)
    y = int(np.clip(ys.mean() - rng.uniform(24, TILE - 24), 0, image.shape[0] - TILE))
    x = int(np.clip(xs.mean() - rng.uniform(24, TILE - 24), 0, image.shape[1] - TILE))
    image, masks = image[y:y + TILE, x:x + TILE], masks[:, y:y + TILE, x:x + TILE]
    if rng.random() < .5:
        image, masks = image[:, ::-1], masks[:, :, ::-1]
    if rng.random() < .5:
        image, masks = image[::-1], masks[:, ::-1]
    turns = int(rng.integers(4))
    image = np.rot90(image, turns)
    masks = np.rot90(masks, turns, axes=(1, 2))
    return image, masks


def load_checkpoint(path, model):
    checkpoint = torch.load(path, map_location="cpu", weights_only=True)
    if checkpoint.get("schema") != MODEL_SCHEMA:
        raise ValueError("This checkpoint is not an XRI spot Mask R-CNN model.")
    model.load_state_dict(checkpoint["model"])
    return checkpoint["training"]


def train(args):
    if args.synthetic < 0 or args.fine_tune < 0 or args.batch < 1 or args.synthetic + args.fine_tune < 1:
        raise ValueError("Supply positive training counts and batch size.")
    if args.synthetic == 0 and not args.resume:
        raise ValueError("Fine tuning requires a pretrained checkpoint.")
    rng = np.random.default_rng(args.seed)
    torch.manual_seed(args.seed)
    device = torch.device(args.device)
    model = build_model().to(device)
    history = load_checkpoint(args.resume, model) if args.resume else []
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    def run(count, phase, sample, rate):
        if not count:
            return
        model.train()
        optimizer = torch.optim.AdamW(model.parameters(), lr=rate)
        running = 0.
        for first in range(0, count, args.batch):
            batch = [training_sample(*sample()) for _ in range(min(args.batch, count - first))]
            images = [image.to(device) for image, _ in batch]
            targets = [{key: value.to(device) for key, value in target.items()} for _, target in batch]
            loss = sum(model(images, targets).values())
            if not torch.isfinite(loss):
                raise ValueError("Training encountered a nonfinite loss.")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.)
            optimizer.step()
            running += loss.item()
            step = first // args.batch + 1
            if step % 25 == 0 or first + len(batch) == count:
                print(f"{phase} {first + len(batch)}/{count} mean loss {running / step:.4f}", flush=True)
            if (first + len(batch) == count) or step % 100 == 0:
                temp = output.with_suffix(output.suffix + ".tmp")
                torch.save({"schema": MODEL_SCHEMA, "model": model.state_dict(),
                            "training": history + [{"phase": phase, "images": first + len(batch),
                                                    "mean_loss": running / step, "seed": args.seed}]}, temp)
                temp.replace(output)
        history.append({"phase": phase, "images": count, "mean_loss": running / step, "seed": args.seed})

    run(args.synthetic, "synthetic", lambda: synthetic(rng), .001)
    if args.fine_tune:
        examples = list(real_crops(args))
        if not examples:
            raise ValueError("Fine tuning requires real masks or valid ROI boxes and a real detector scan.")
        phase = "real_masks_and_weak_boxes" if args.real_masks and args.real_scan else (
            "real_masks" if args.real_masks else "real_weak_boxes")
        run(args.fine_tune, phase, lambda: crop_real(rng, examples), .0002)


def positions(size):
    last = max(0, size - TILE)
    return sorted(set(range(0, last + 1, TILE - 56)) | {last})


def detect(model, raw, score, mask_threshold):
    height, width = raw.shape
    candidates = []
    image = features(raw)
    with torch.inference_mode():
        for y in positions(height):
            for x in positions(width):
                tile = image[:, y:y + TILE, x:x + TILE]
                _, h, w = tile.shape
                padded = np.pad(tile, ((0, 0), (0, TILE - h), (0, TILE - w)))
                result = model([torch.from_numpy(padded)])[0]
                for confidence, probability in zip(result["scores"].numpy(), result["masks"][:, 0].numpy()):
                    if confidence < score:
                        continue
                    mask = probability[:h, :w] >= mask_threshold
                    labels, count = ndimage.label(mask)
                    if not count:
                        continue
                    sizes = np.bincount(labels.ravel())
                    sizes[0] = 0
                    ys, xs = np.nonzero(labels == sizes.argmax())
                    if xs.size < 8:
                        continue
                    # Prefer complete masks where neighboring tiles overlap.
                    clipped = ((x > 0 and xs.min() == 0) or (y > 0 and ys.min() == 0) or
                               (x + w < width and xs.max() == w - 1) or
                               (y + h < height and ys.max() == h - 1))
                    candidates.append((bool(clipped), float(confidence), (ys + y) * width + xs + x))
    selected = []
    buckets = {}
    for _, confidence, pixels in sorted(candidates, key=lambda item: (item[0], -item[1])):
        ys, xs = np.divmod(pixels, width)
        cells = [(y, x) for y in range(ys.min() // 64, ys.max() // 64 + 1)
                 for x in range(xs.min() // 64, xs.max() // 64 + 1)]
        nearby = {index for cell in cells for index in buckets.get(cell, ())}
        if any(np.intersect1d(pixels, selected[index][1], assume_unique=True).size
               >= .55 * min(pixels.size, selected[index][1].size) for index in nearby):
            continue
        for cell in cells:
            buckets.setdefault(cell, []).append(len(selected))
        selected.append((confidence, pixels))
    return selected


def measure(raw, confidence, pixels):
    height, width = raw.shape
    y, x = np.divmod(pixels, width)
    bounds = [int(x.min()), int(y.min()), int(x.max() + 1), int(y.max() + 1)]
    x1, y1, x2, y2 = bounds
    left, top, right, bottom = max(0, x1 - 8), max(0, y1 - 8), min(width, x2 + 8), min(height, y2 + 8)
    mask = np.zeros((bottom - top, right - left), dtype=bool)
    mask[y - top, x - left] = True
    ring = ndimage.binary_dilation(mask, iterations=7) & ~ndimage.binary_dilation(mask, iterations=2)
    ring_values = raw[top:bottom, left:right][ring]
    finite_ring = ring_values[np.isfinite(ring_values)]
    background = float(np.median(finite_ring)) if finite_ring.size else 0.
    values = np.where(np.isfinite(raw[y, x]), raw[y, x], background)
    weights = np.maximum(values - background, 0)
    xy = np.column_stack((x, y)).astype(np.float64)
    center = np.average(xy, axis=0, weights=weights) if weights.sum() > 0 else xy.mean(axis=0)
    geometry = xy - xy.mean(axis=0)
    eigenvalues, axes = np.linalg.eigh(geometry.T @ geometry / len(x))
    delta = np.diff(pixels, prepend=-2, append=int(pixels[-1]) + 2)
    starts = np.flatnonzero(delta[:-1] != 1)
    ends = np.flatnonzero(delta[1:] != 1) + 1
    return {
        "score": confidence, "bbox_xyxy": bounds, "mask_area": int(x.size),
        "centroid_xy": center.tolist(), "integrated_counts": float(values.sum()),
        "net_counts": float((values - background).sum()), "background_per_pixel": background,
        "ellipse_major_px": float(4 * np.sqrt(max(eigenvalues[1], 0))),
        "ellipse_minor_px": float(4 * np.sqrt(max(eigenvalues[0], 0))),
        "ellipse_angle_deg": float(np.degrees(np.arctan2(axes[1, 1], axes[0, 1])) % 180),
        "mask_runs": [[int(pixels[i]), int(j - i)] for i, j in zip(starts, ends)],
    }


def infer(args):
    if not 0 < args.score <= 1 or not 0 < args.mask_threshold <= 1 or args.frame < -1:
        raise ValueError("Thresholds must be in (0,1], with frame index at least -1.")
    model = build_model()
    history = load_checkpoint(args.weights, model)
    model.eval()
    print("Reading detector data", flush=True)
    raw = read_detector(args.file, args.detector, args.frame)
    if not np.isfinite(raw).any():
        raise ValueError("The selected detector data has no finite pixels.")
    print("Running Mask R-CNN", flush=True)
    candidates = detect(model, raw, args.score, args.mask_threshold)
    measured = [(measure(raw, confidence, pixels), pixels) for confidence, pixels in candidates]
    measured.sort(key=lambda item: (item[0]["centroid_xy"][1], item[0]["centroid_xy"][0]))
    gray = np.rint(features(raw)[2] * 200).astype(np.uint8)
    overlay = np.repeat(gray[..., None], 3, axis=2)
    colors = [(34, 219, 192), (242, 185, 60), (233, 108, 176), (113, 188, 255)]
    for index, (spot, pixels) in enumerate(measured, 1):
        spot["id"] = index
        color = np.array(colors[(index - 1) % len(colors)])
        overlay.reshape(-1, 3)[pixels] = .45 * overlay.reshape(-1, 3)[pixels] + .55 * color
    preview = Image.fromarray(overlay)
    draw = ImageDraw.Draw(preview)
    for spot, _ in measured:
        x, y = spot["centroid_xy"]
        draw.text((x + 3, y + 2), str(spot["id"]), fill="white", stroke_width=1, stroke_fill="black")
    preview.save(args.preview)
    report = {"schema": "xri.roi.maskrcnn.result.v1", "source": str(Path(args.file).resolve()),
              "detector": args.detector, "frame": args.frame, "shape_yx": list(raw.shape),
              "score_threshold": args.score, "mask_threshold": args.mask_threshold, "training": history,
              "mask_encoding": "zero-based row-major runs [start,length]",
              "spots": [spot for spot, _ in measured]}
    Path(args.output).write_text(json.dumps(report, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    print(f"Detected {len(measured)} spots", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=4)
    commands = parser.add_subparsers(dest="command", required=True)
    training = commands.add_parser("train")
    training.add_argument("--synthetic", type=int, default=3000)
    training.add_argument("--fine-tune", type=int, default=0)
    training.add_argument("--real-masks", nargs="*", help="NPZ files containing image[H,W] and masks[N,H,W]")
    training.add_argument("--real-scan")
    training.add_argument("--roi-boxes", help="XRI ROI JSON for generating weak mask labels")
    training.add_argument("--resume")
    training.add_argument("--batch", type=int, default=2)
    training.add_argument("--seed", type=int, default=1219)
    training.add_argument("--device", default="cpu")
    training.add_argument("--output", required=True)
    inference = commands.add_parser("infer")
    inference.add_argument("--file", required=True)
    inference.add_argument("--detector", required=True, choices=DETECTORS)
    inference.add_argument("--frame", type=int, default=-1)
    inference.add_argument("--weights", required=True)
    inference.add_argument("--score", type=float, default=.55)
    inference.add_argument("--mask-threshold", type=float, default=.5)
    inference.add_argument("--output", required=True)
    inference.add_argument("--preview", required=True)
    args = parser.parse_args()
    torch.set_num_threads(max(1, args.threads))
    torch.set_num_interop_threads(1)
    (train if args.command == "train" else infer)(args)


if __name__ == "__main__":
    main()

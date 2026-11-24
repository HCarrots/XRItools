"""Export the optical calibration needed by XRI from a completed native-band cache."""

import argparse
import hashlib
import json
import tempfile
from pathlib import Path

import h5py


def export_preset(cache: Path, config_path: Path, output: Path) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    source = Path(config["input"]["calibration_nxs"])
    roi = Path(config["input"]["roi_table"])
    if source.name != "target_kapton_55.nxs" or config["spatial_calibration"]["match_sample_band"]:
        raise ValueError("Use the target_kapton_55 native-band calibration configuration.")

    output.mkdir(parents=True, exist_ok=True)
    with (
        h5py.File(cache / "calibration.h5", "r") as raw,
        h5py.File(cache / "psf.h5", "r") as models,
        tempfile.TemporaryDirectory(dir=output) as temporary,
    ):
        if Path(raw.attrs["source_file"]).resolve() != source.resolve() or Path(raw.attrs["roi_file"]).resolve() != roi.resolve():
            raise ValueError("The cache references a different foil scan or ROI scheme.")
        acquisition = config["acquisition"]
        for field, attr in (("monitor", "monitor"), ("calibration_exposure_s", "exposure_s"),
                            ("foil_thickness_mm", "foil_thickness_mm")):
            if acquisition[field] != raw.attrs[attr]:
                raise ValueError(f"The cache uses a different {field}.")
        roi_scheme = json.loads(roi.read_text(encoding="utf-8-sig"))
        roi_scheme.pop("legacy_source_files", None)
        roi_scheme.pop("legacy_conversion_notes", None)
        roi_scheme["source"] = source.name
        manifest = json.loads(models["manifest"][()])
        if not manifest["models"]:
            raise ValueError("The cache has no fitted optical models.")
        with source.open("rb") as scan:
            digest = hashlib.file_digest(scan, "sha256").hexdigest()
        preset = {
            "schema": "XRI.prepared-foil.v1",
            "source_scan": source.name,
            "source_sha256": digest,
            "detector": acquisition["detector"],
            "monitor": acquisition["monitor"],
            "calibration_exposure_s": acquisition["calibration_exposure_s"],
            "foil_thickness_mm": acquisition["foil_thickness_mm"],
            "detector_pixel_pitch_mm": acquisition["detector_pixel_pitch_mm"],
            "psf": config["psf"],
            "loop_order": ["X", "DCM"],
            "actual_frames": len(raw["frames/i0"]),
        }
        _, height, width = raw[f"detectors/{preset['detector']}/band_response_xrc"].shape
        preset["image_shape"] = [height, width]
        staged = Path(temporary)
        name = "target_kapton_55"
        data_path = staged / f"{name}.calibration.h5"
        roi_path = staged / f"{name}.roi.json"
        with h5py.File(data_path, "w", libver="earliest") as result:
            for key, value in raw.attrs.items():
                if key not in {"source_file", "roi_file"}:
                    result.attrs[key] = value
            result.attrs["source_file"] = source.name
            result.attrs["roi_file"] = roi_path.name
            result.create_dataset("xri_prepared_foil", data=json.dumps(preset, separators=(",", ":")))
            del manifest["key"]
            result.create_dataset("manifest", data=json.dumps(manifest, separators=(",", ":")))
            for key in ("x_mm", "energy_kev", "energy_weights_kev"):
                raw.copy(key, result)
            for key in sorted(manifest["models"].keys() | manifest["rejected"].keys()):
                raw_group = raw[f"rois/{key}"]
                target = result.create_group(f"rois/{key}")
                for attr, value in raw_group.attrs.items():
                    target.attrs[attr] = value
                for item in ("pixels_rc", "bbox_rc", "R_band_px", "R_band_variance_px"):
                    raw_group.copy(item, target)
            if "templates" in models:
                models.copy("templates", result)
        roi_path.write_text(json.dumps(roi_scheme, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        data_path.replace(output / data_path.name)
        roi_path.replace(output / roi_path.name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    export_preset(args.cache, args.configuration, args.output)

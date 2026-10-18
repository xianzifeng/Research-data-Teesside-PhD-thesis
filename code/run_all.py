"""Run the full Python replication pipeline in order and tee all console output to outputs/run_log.txt."""
import subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "outputs"; OUT.mkdir(exist_ok=True)
STEPS = ["00_generate_synthetic_data.py", "01_data_preprocessing.py", "02_measurement_checks.py",
         "03_structural_paths_crosscheck.py", "04_bootstrap_mediation.py", "05_iv_estimation.py",
         "06_moderation_interaction.py", "07_robustness_checks.py", "08_mmdv_worksheet.py",
         "09_prepare_mplus_data.py"]
# Skip synthetic generation when the (undistributed) real panel is present
if (HERE.parent / "data" / "real" / "huafeng_panel_person_wave.csv").exists():
    STEPS.remove("00_generate_synthetic_data.py")

with open(OUT / "run_log.txt", "w", encoding="utf8") as log:
    for s in STEPS:
        t0 = time.time()
        print(f"\n>>> {s}")
        r = subprocess.run([sys.executable, str(HERE / s)], capture_output=True, text=True, cwd=HERE)
        print(r.stdout); log.write(f"\n>>> {s}\n{r.stdout}")
        if r.returncode:
            print(r.stderr); log.write(r.stderr); sys.exit(f"{s} failed")
        print(f"    ({time.time() - t0:.1f} s)")
print(f"\nAll steps completed. Log: {OUT / 'run_log.txt'}")

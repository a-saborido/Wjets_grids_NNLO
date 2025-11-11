import os
import numpy as np
import subprocess
from pathlib import Path

# Define directory containing .tab.gz files
THISDIR = os.environ.get("THISDIR", ".")  # Default to current directory if not set

# Generate 100 scale variation values between 0.1 and 32
#scale_variations = np.linspace(0.1, 10, 100)
scale_variations = np.logspace(-1.5, 1.5, 100)

# Get list of .tab.gz files
tab_files = list(Path(THISDIR).glob("NNLO*.tab.gz"))

for file_path in tab_files:
    for scale in scale_variations:
        # Prepare output filename
        out_filename = f"{file_path}_{scale:.5f}dscale1.txt.out"
        out_path = Path(out_filename)

        # Touch the output file (create it if it doesn't exist)
        out_path.touch()

        # Build the command
        cmd = [
            "fnlo-tk-cppread",
            str(file_path),
            "NNPDF31_nnlo_as_0118",
            "3",
            "LHAPDF",
            "no",
            "kScale1",
            f"{scale}"
        ]

        # Run the command and write output to file
        with open(out_path, "w") as outfile:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.STDOUT)

        print(f"finished processing {out_path}")

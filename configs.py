import os
import argparse

# ========== General configs ==========
current_dir = os.path.dirname(__file__)
data_path = os.path.join(current_dir, "data")
results_path = os.path.join(current_dir, "results")

# ========== Args Parser ==============
parser = argparse.ArgumentParser()

parser.add_argument(
    '-m', '--model',
    type=str,
    default="nearest_neighbors",
    help="Model name"
)

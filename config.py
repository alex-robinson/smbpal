"""
Script to generate Makefile with the correct compiler
configuration for a given machine.
"""

from subprocess import *
import argparse
import sys
import re

def parse_makefile(filename="Makefile"):
    """
    Parses a Makefile and prints a list of available `make` commands (targets).
    
    Args:
        filename (str): The path to the Makefile. Defaults to "Makefile".
    """
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            
        targets = []
        for line in lines:
            # Strip comments and whitespace
            line = line.split('#')[0].strip()
            
            # Match lines with targets (e.g., target: dependencies)
            match = re.match(r'^([a-zA-Z0-9_-]+)\s*:\s*(.*)', line)
            if match and not line.startswith('\t'):
                target = match.group(1)
                targets.append(target)
        
        print("Available `make` commands (targets):")
        for target in targets:
            print(f"  - {target}")
        print("")
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Manage command-line arguments to load compiler option
parser = argparse.ArgumentParser()
parser.add_argument(
    "config",
    metavar="CONFIG",
    type=str,
    help="Name of config file to use (e.g., config/aire_gfortran)",
)
args = parser.parse_args()


# Determine locations
target_dir = "./"
config_dir = "config/"
config_path = args.config

# Load template Makefile and insert compiler configuration
makefile0 = open(config_dir + "Makefile").read()
compile_info = open(config_path).read()
makefile1 = makefile0.replace("<COMPILER_CONFIGURATION>", compile_info)

# Write the new Makefile to the main directory
open(target_dir + "Makefile", "w").write(makefile1)

print(
    "".join(
        [
            "\n\nMakefile configuration complete for configuration file: ",
            config_path,
            "\n",
        ]
    )
)

# Print usage
parse_makefile("Makefile")

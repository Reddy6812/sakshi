import os
import shutil
import random

def move_random_n_files(source_dir, destination_dir, n):
    # Check if source and destination directories exist
    if not os.path.exists(source_dir):
        print(f"Source directory {source_dir} does not exist.")
        return
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)
        print(f"Destination directory {destination_dir} did not exist. Created it.")

    # Get a list of files in the source directory
    files = [f for f in os.listdir(source_dir) if os.path.isfile(os.path.join(source_dir, f))]

    # Check if there are enough files to move
    if len(files) < n:
        print(f"Not enough files to move. Found only {len(files)} files.")
        n = len(files)  # Adjust n to the number of available files

    # Randomly select n files from the list
    random_files = random.sample(files, n)

    # move randomly selected files
    for file_name in random_files:
        file_path = os.path.join(source_dir, file_name)
        shutil.move(file_path, destination_dir)
        print(f"Moved {file_name} to {destination_dir}")


source_directory = '/blue/yanjun.li/vi.gade1/seabra-li/data/usefortrain/'
destination_directory = '/blue/yanjun.li/vi.gade1/seabra-li/data/enaminesdf2pdbqt/'
n_files = 2000  # Number of files to move

move_random_n_files(source_directory, destination_directory, n_files)

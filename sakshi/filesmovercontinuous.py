import os
import shutil

def move_first_n_files(source_dir, destination_dir, n):
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

    # move first n files
    for i in range(n):
        file_path = os.path.join(source_dir, files[i])
        shutil.move(file_path, destination_dir)
        print(f"Moved {files[i]} to {destination_dir}")


source_directory = '/blue/yanjun.li/vi.gade1/seabra-li/data/enaminesdf2pdbqt/'
destination_directory = '/blue/yanjun.li/vi.gade1/seabra-li/data/enaminemoved/'
n_files = 400000  # Number of files to move

move_first_n_files(source_directory, destination_directory, n_files)

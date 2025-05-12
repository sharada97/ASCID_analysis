import os
import glob

def get_folder_size(folder_path):
    """Calculate the total size of all files in a folder."""
    total_size = 0
    for dirpath, _, filenames in os.walk(folder_path):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            # Skip if the file is broken or inaccessible
            if os.path.exists(filepath):
                total_size += os.path.getsize(filepath)
    return total_size

def filter_folders_by_size(folder_list, size_threshold):
    """
    Filter folders by their total memory usage.

    Args:
        folder_list (list): List of folder paths.
        size_threshold (int): Minimum folder size in bytes.

    Returns:
        list: List of folder paths meeting the size criteria.
    """
    large_folders = [folder for folder in folder_list if get_folder_size(folder) > size_threshold]
    return large_folders

# Example usage
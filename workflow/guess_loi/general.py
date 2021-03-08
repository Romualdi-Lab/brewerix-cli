from os import remove


def remove_files(files):
    if isinstance(files, str):
        safe_remove(files)
    else:
        for file in files:
            safe_remove(file)


def safe_remove(file):
    try:
        remove(file)
    except OSError:
        pass


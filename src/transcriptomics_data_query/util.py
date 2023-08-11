import os.path
import tarfile

def is_valid_tar_member(member: tarfile.TarInfo, target_dir: str) -> bool:
    """Check if a tar member is safe to extract. `target_dir` should be an absolute path."""
    member_dest = os.path.abspath(os.path.join(target_dir, member.name))
    return member_dest.startswith(target_dir)

def extract_tar(tar_file, target_dir, delete_tar=False):
    """
    Extract a tar file to a target directory.

    Parameters
    ----------
    tar_file : str
        Path to the tar file.
    target_dir : str
        Path to the target directory.
    delete_tar : bool, optional
        If True, delete the tar file after extraction, default is False.

    """
    target_dir = os.path.abspath(target_dir)

    with tarfile.open(tar_file, 'r') as tar:
        safe_members = [member for member in tar.getmembers()
                        if is_valid_tar_member(member, target_dir)]
        tar.extractall(path=target_dir, members=safe_members)

    if delete_tar:
        os.remove(tar_file)

import os
from pathlib import Path
import argparse

def update_ncbi_email(new_email, pkg_path):
    """
    Update the email in the designated text file.

    Parameters
    ----------
    new_email : str
        The new email address to be written to the file.
    package_path : str
        The path to the directory containing the email text file.

    """
    email_file_path = os.path.join(pkg_path, 'email_for_ncbi_tracking.txt')
    with open(email_file_path, 'w', encoding='utf-8') as file:
        file.write(new_email)
    print(f"Email updated to {new_email}")

def main():
    parser = argparse.ArgumentParser(description="Update email for NCBI tracking.")
    parser.add_argument("email", help="New email address to set.")
    args = parser.parse_args()

    try:
        package_path = Path(__file__).resolve().parent
    except NameError:
        # Jupyter Notebook
        from IPython import get_ipython
        package_path = Path(get_ipython().getoutput('pwd')[0])

    update_ncbi_email(args.email, package_path)

if __name__ == "__main__":
    main()

import os
import re

def install_packages_from_md(file_path):
    with open(file_path, "r") as file:
        content = file.read()

    # Extract package names and versions from the Markdown content
    package_pattern = r"\s*-\s*(\S+)==(\S+)\s*"
    packages = re.findall(package_pattern, content)

    # Install each package using pip
    for package, version in packages:
        install_command = f"pip install {package}=={version}"
        os.system(install_command)

if __name__ == "__main__":
    md_file_path = "python_packages.md"
    install_packages_from_md(md_file_path)

import os
import pkg_resources


def get_python_packages():
    # Get a list of all installed Python packages and their versions
    installed_packages = [f"{pkg.key}=={pkg.version}" for pkg in pkg_resources.working_set]
    installed_packages.sort()
    return installed_packages


def create_markdown_file(file_path, python_version, packages):
    # Create the Markdown content
    markdown_content = f"# Python Packages and Version\n\n"
    markdown_content += f"Python Version: {python_version}\n\n"
    markdown_content += "Installed Packages:\n\n"

    for package in packages:
        markdown_content += f"- {package}\n"

    # Write the Markdown content to the file
    with open(file_path, "w") as file:
        file.write(markdown_content)


if __name__ == "__main__":
    # Get Python version
    python_version = f"{os.sys.version_info.major}.{os.sys.version_info.minor}.{os.sys.version_info.micro}"

    # Get installed Python packages
    packages = get_python_packages()

    # Create the Markdown file
    create_markdown_file("python_packages.md", python_version, packages)

    print("Python packages and version information have been written to 'python_packages.md'")


import subprocess

def get_gams_version():
    try:
        # Run the 'gams --version' command and capture the output
        result = subprocess.run(['gams', '--version'], capture_output=True, text=True, check=True)

        # Extract the version information from the output
        version = result.stdout.strip()
        return version
    except FileNotFoundError:
        return "GAMS is not installed or not added to PATH."
    except subprocess.CalledProcessError as e:
        return f"Error while checking GAMS version: {e}"

if __name__ == "__main__":
    gams_version = get_gams_version()
    print("GAMS Version:", gams_version)

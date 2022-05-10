import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="macrodna-xiru.h", # Replace with your own username
    version="0.0.1",
    author="Mohammadamin Edrisi, Xiru Huang",
    author_email="xiru.huang@rice.edu",
    description="A high-accuracy tool for integrating scDNA-seq and scRNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xiru-huang/MaCroDNA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
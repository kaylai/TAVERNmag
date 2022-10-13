import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="TAVERNmag",
    version="0.2.0",
    author="Kayla Iacovino",
    author_email="kaylaiacovino@gmail.com",
    description=("A thermodynamic model code for magmatic volatiles."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaylai/TAVERNmag",
    packages=setuptools.find_packages(),
    install_requires=[
            'pandas',
            'numpy',
            'scipy',
            'sympy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)

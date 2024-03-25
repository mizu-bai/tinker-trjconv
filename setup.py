import setuptools

setuptools.setup(
    name="tinker-trjconv",
    version="0.1.0",
    author="mizu-bai",
    author_email="shiragawa4519@outlook.com",
    description="A Tinker trajectory converter",
    url="https://github.com/mizu-bai/tinker-trjconv",
    packages=setuptools.find_packages(),
    classifiers=[],
    python_requires=">=3.8",
    install_requires=[
        "numpy>1.20.0",
    ],
)

# DAESIM2 Analysis

## Description

This Python package supports the development and testing of the Dynamic Agro-Ecosystem Simulator 2 (DAESIM2), a biophysical model for agro-ecosystems. The package provides tools for conducting sensitivity analysis, calibration, and uncertainty analysis, helping to evaluate and improve model performance.

References: Taghikhah et al. (2022) https://doi.org/10.1016/j.ecolmodel.2022.109930

## Installation

TODO

### Installing tools on Mac

TODO

## Getting started

Make your conda environment called `daesim2-analysis` using the environment.yml file and the command:

```sh
conda env create --name daesim2-analysis --file environment.yml
```

Then, you can activate your new environment with the command:

```sh
conda activate daesim2-analysis
```

Now you have created the conda virtual environment to work in, you must install the package. To install the package use the command:

```sh
pip install -e .
```

Hopefully everything installed successfully. 

You can now run the code or open a jupyter notebook to start testing. To open a jupyter notebook run `jupyter-notebook` to open a Jupyter Notebook server or run `jupyter-lab` to open a JupyterLab server (recommended). From there, go into the notebooks directory and work through the examples.

Once you're finished, you can deactivate the conda environment with the command:

```sh
conda deactivate
```

### Linking to the DAESIM2 GitHub repository

This repository is linked to the DAESIM2 biophysical model source code [DAESIM2](https://github.com/NortonAlex/DAESIM). The analysis tools are kept separate from the biophysical model code for better organization and maintainability. Keeping them independent prevents unnecessary complexity in the model codebase and ensures that the analysis tools can evolve separately without impacting core model functionality.

This separation benefits both model users and developers:

* Users can run and interact with the DAESIM2 model without needing to understand or install additional analysis tools (including their specialized Python libraries).

* Developers can work on analysis tasks (e.g., calibration, sensitivity analysis, uncertainty analysis) without interfering with the core model implementation, also allowing for a stable reference to DAESIM2 model versions.


#### Updating the Environment and Dependencies

If the DAESIM2 biophysical model (linked via GitHub) has been updated or if you want to reference a difference version of the model, the following steps must be followed to ensure those changes are properly integrated into `daesim2-analysis`:

	1. Get DAESIM2 commit hash

	Note the GitHub commit hash of the DAESIM2 version that you want to integrate 

	2. Update the dependencies in daesim2-analysis

		2.1 Update environment.yml

		Locate the - pip: section in environment.yml and update the git URL reference with the new commit hash.

  		2.2 Update setup.cfg

		Locate the install_requires section within the setup.cfg file and update the git URL reference with the new commit hash. 

	3. Update the Conda Environment

		Option 1: Update the Existing Environment (Recommended)

		Run the following command:

		```sh
		conda env update --file environment.yml --prune
		```

		This will update the environment while keeping existing dependencies intact.

		Option 2: Recreate the Environment (If Issues Occur)

		Sometimes, conda env update does not reflect the new changes and the environment must be removed and created again:

		```sh
		conda env remove -n daesim2-analysis
		conda env create -f environment.yml
		conda activate daesim2-analysis
		```

	4. Reinstall the Package

		After updating dependencies, reinstall `daesim2-analysis` to ensure the latest changes are recognized:

		```sh
		pip install -e .
		```

		To force the latest GitHub dependency installation:

		```sh
		pip install --force-reinstall git+https://github.com/NortonAlex/DAESIM@<new-commit-hash>#egg=daesim
		```

	5. Verify the Update

		To ensure everything is correctly installed:

		Open a Python shell or Jupyter Notebook.

		Try importing the packages:

		```sh
		import daesim
		import daesim2_analysis
		```

		If there are issues, verify the installed package versions:

		```sh
		pip list | grep daesim
		```

		If any errors occur, consider clearing the pip cache and reinstalling:

		```sh
		pip cache purge
		pip install -e .
		```

## General guidance on working in this repository

### Repository guiding principles

Text todo: Outline how to use this repository. e.g. /data/ should only contain a test set of forcing data, do not commit large datasets to it. Similarly, do not commit anything to /results/. Keep the repo clean and tidy. 

### Branch etiquette

In general, don't push to other people's branches (things can get weird if two people work on the same branch at once).
Instead, check out the branch and then create your own branch in case you want to make any commits (see "Checking out a branch related to a merge request" below).
Then, to add your changes, make a merge request from your own branch back into the original branch.

## Support

If you need support, the first place to go is the [issue tracker] (todo: link the repository / issues).
From there, you can tag other model users to ask for help.
As a second step, reach out directly to your collaborators.

## Other helpful background

### Why do we use Anaconda?

There are benefits to using Anaconda rather than just using pip. As discussed in this article (https://pythonspeed.com/articles/conda-vs-pip/) the main benefits are portability and reproducibility.
- Portability across operating systems: Instead of installing Python in three different ways on Linux, macOS, and Windows, you can use the same environment.yml on all three.
- Reproducibility: It’s possible to pin almost the whole stack, from the Python interpreter upwards.
- Consistent configuration: You don’t need to install system packages and Python packages in two different ways; (almost) everything can go in one file, the environment.yml.

Using Conda also addresses another problem: How to deal with Python libraries that require compiled code? At some point we may like to convert some of the DAESIM2 source code into another, more computationally efficient language (e.g. C++ or Fortran) so that we can run larger simulations. If we decide to do that we will need to compile the source code in a language other than Python, in which case pip would be insufficient and Conda will be required. 

### Why do I only see `.py` files in the `notebooks` directory?

Tracking notebooks is a pain because of all the extra metadata that is in them.
As a result, we use [jupytext](https://jupytext.readthedocs.io/).
This allows us to track clean `.py` files while keeping all the power of jupyter notebooks.
You'll notice that we ignore `*.ipynb` in our `.gitignore`.

Thanks to jupytext, you can just click on the `.py` files and they will open like standard notebooks.

### Checking out a branch related to a merge request

```sh
# Fetch any changes on the remote server
git fetch origin
# If you get no output, you're up to date
# If you get output, it will show you which branches have changed
# on the remote server

# Checkout the others' branch and create a local branch to track it
# git checkout -b branch-name origin/branch-name
# for example:
git checkout -b local-source origin/local-source

# Activate your environment
conda activate daesim2-analysis

[remove]# Call make (just in case, you don't always have to do this)
# make conda-environment

# If any dependencies have changed, update your conda environment
conda env update --file environment.yml --prune

# Checkout (create) your own branch in case you want to make any commits
# (typically we use the branch name plus our initials)
git checkout -b test-notebook-an

# You're now ready to work
# E.g. by starting a notebook server and looking at the notebooks
jupyter notebook
```

## Reading and viewing

Not all of these will be neccesary for everyone, but a helpful list

### What is self

https://www.educative.io/answers/what-is-self-in-python

### General approach to coding

- [Clean Code](https://thixalongmy.haugiang.gov.vn/media/1175/clean_code.pdf) (buying the book is also a good option)
- [Refactoring](http://silab.fon.bg.ac.rs/wp-content/uploads/2016/10/Refactoring-Improving-the-Design-of-Existing-Code-Addison-Wesley-Professional-1999.pdf) (buying the book is also a good option)
- [Refactoring guru](refactoring.guru), incredible resource for understanding coding patterns and how to make them better (There is also a book, could be worth investing in)
- [End of object inheritance](https://www.youtube.com/watch?v=3MNVP9-hglc) This one is hard to explain and understand until you start writing lots of code, but its worth watching (and re-watching) to understand the coding style you see in new climate models (e.g. MAGICC)
- [End of object inheritance](https://www.youtube.com/watch?v=3MNVP9-hglc) This one is hard to explain and understand until you start writing lots of code, but its worth watching (and re-watching) to understand the coding style you see in new climate models (e.g. MAGICC)
- [Composition over inheritance](https://www.youtube.com/watch?v=0mcP8ZpUR38) A nice explainer to see the principles of the above in practice
- [Dependency injection vs. inversion](https://www.youtube.com/watch?v=2ejbLVkCndI) A further explainer to see the above in practice

### Scientific software

- [Research software engineering with Python](https://merely-useful.tech/py-rse/) A little bit out of date, but a good resource for general practice and examples of developing software with python

### Numerical coding

- [Numerical recipes](http://numerical.recipes/book/book.html) (buying the book is also a good option)

### Big climate model output and handling netCDF files

- [Software Carpentry's introduction to Python for atmosphere and ocean scientists](https://carpentries-lab.github.io/python-aos-lesson/), particularly if you think you're going to be working with netCDF files
- [CMIP6 in the cloud](https://medium.com/pangeo/cmip6-in-the-cloud-five-ways-96b177abe396), particularly if you're going to be dealing with CMIP data a lot (although speak with someone about how far you need to go, it's unlikely that every step will be relevant)

### Miscellaneous

- [Basic introduction to Jupyter notebooks](https://realpython.com/jupyter-notebook-introduction/)

## Recommended tools

- An IDE (e.g. pycharm) or more light-weight editor (e.g. Sublime)
- If on Mac, homebrew
- Git
- Make
- Slack
- gfortran
- cmake

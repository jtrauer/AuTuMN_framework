# AuTuMN_framework
Contains examples of the general framework of the (private) AuTuMN platform to illustrate modular structure.

Each module that we have coded in the full AuTuMN codebase is represented in this repository and each of the major
functions of the module is demonstrated, including each module represented as a rectangle in Figure 2 of our manuscript
describing the operation of AuTuMN (with the exception of the base module, which we have released in full). Each module
is fully functional such that the codebase runs as a whole using each component module. Therefore, AuTuMN_framework
could be used to undertake relatively simple country-level analyses and could rapidly be developed to undertake
AuTuMN-like analyses.

We believe this codebase will actually be considerably more useful to TB modellers than releasing our full codebase in
its current form. This is because, although AuTuMN is extensively coded and commented, we have not produced user
documentation that would allow programmers previously unfamiliar with our platform to fully interact with it and
to make further modifications. For example, the curve fitting module produces piecewise polynomial spline functions
for use by the data processing module, but has several working approaches to achieve this, not all of which are
currently in use or fully explained to unfamiliar users. Similarly, our current GUI is currently being replaced by a
more advanced version.

The framework repository includes the following modules:
•	Interface module
    o	Presents all the data structures for the user to interact with the remaining codebase and the basis for creating
        more visually appealing graphical user interfaces
•	Model runner module
    o	Presents the basis for running the platform with multiple different purposes, with scenario analysis and
        uncertainty presented as examples (which are purposes that were fully coded by our team and have minimal
        dependencies)
    o	Calculates consistently formatted output structures
•	Data processing module
    o	Formats both fixed and time-variant input parameters for use by the TB model module
•	Spreadsheet module
    o	Reads one external spreadsheet (the WHO’s Global TB Report), making the extracted data available to the data
        processing module
•	Curve fitting module
    o	Fits a sinusoidal piecewise function to input data
•	Outputs module
    o	Produces graphical outputs from the data structures created by the model runner
•	Tool kit
    o	Stores static functions for use by any other module
•	TB model
    o	Inherits the BaseModel class from the popdynamics module to create a TB model with one tier of stratification

# Requirements

- python 2.7
- numpy
- scipy
- matplotlib
- platform
- xlrd
- os
- graphviz.py
- graphviz binary on path: http://www.graphviz.org/
- popdynamics binary on path: https://github.com/popdynamics/popdynamics

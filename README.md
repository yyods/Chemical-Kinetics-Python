# Chemical Kinetics with Python and SymPy

This repository contains Python code and resources accompanying the teaching material for _Introduction to Python and SymPy for Chemical Kinetics_. The purpose of this repository is to help chemistry students follow along with the teaching material, understand Python and SymPy, and utilize Jupyter Notebook for solving mathematical problems related to chemical kinetics.

## 📚 About This Repository

The repository provides prepared Python scripts that correspond to various sections of the teaching material. These scripts help students:

- Automate repetitive calculations and equation solving.
- Symbolically manipulate and analyze kinetic models efficiently.
- Visualize reaction dynamics and interpret chemical behavior using computational tools.

However, all coding exercises in **Section: Hands-On Activities** are not included in this repository, as students are encouraged to practice and develop their own solutions to build familiarity with coding.

## 📚 Repository Structure

The repository is organized into separate folders corresponding to the chapters of the teaching material:

```
📆 Chemical-Kinetics-Python
 ├📂 Chapter1
 ┃ ├📝 overview_python.py
 ┃ ├📝 introduction_sympy.py
 ┃ ├📝 chemical_kinetics_review.py
 ┃ ├📝 jupyter_notebook_basics.py
 ┃ └📝 sympy_functions_tutorial.py
 ├📂 Chapter2
 ├📂 Chapter3
 ├...
 └📝 README.md
```

- Each chapter folder contains Python scripts related to that chapter.
- All code is designed to be run in Jupyter Notebook for interactive learning.

## 🚀 Getting Started

### 1️⃣ Install Dependencies

Ensure you have Python and the required libraries installed:

```bash
pip install sympy numpy matplotlib scipy jupyter
```

### 2️⃣ Clone the Repository

To download the repository, run:

```bash
git clone https://github.com/yyods/Chemical-Kinetics-Python.git
cd Chemical-Kinetics-Python
```

### 3️⃣ Running Jupyter Notebook

To open Jupyter Notebook and start experimenting:

```bash
jupyter notebook
```

Then, navigate to the relevant `.py` files or convert them into Jupyter Notebook cells for execution.

## 📌 How to Use the Code

Each script in this repository corresponds to a section in the teaching material. You can run the provided scripts, modify them, and explore different parameter values to see how chemical kinetics problems can be modeled computationally.

### Example: Solving a Simple Rate Equation

```python
from sympy import symbols, Eq, solve

# Define symbols
k, A = symbols('k A')

# Define a simple rate equation
rate_eq = Eq(A - k*A, 0)

# Solve for A
solution = solve(rate_eq, A)
print(solution)
```

## 📝 Notes for Students

- **Practice the Hands-On Activities**: The provided code serves as a guide, but it is recommended that students attempt the hands-on exercises themselves.
- **Modify the Code**: Try changing parameter values, adding plots, and experimenting with different kinetics models.
- **Ask Questions**: If you have any questions, feel free to raise an issue on GitHub or discuss in class.

## 🐝 License

This repository is for educational purposes. Feel free to use, modify, and contribute to improve learning.

---

Happy coding! 🚀

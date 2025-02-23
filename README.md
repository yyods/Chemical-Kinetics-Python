# Chemical Kinetics with Python and SymPy

This repository contains Python code and resources accompanying the teaching material for _Introduction to Python and SymPy for Chemical Kinetics_. The purpose of this repository is to help chemistry students follow along with the teaching material, understand Python and SymPy, and utilize Jupyter Notebook for solving mathematical problems related to chemical kinetics.

## 📚 About This Repository

The repository provides prepared Python scripts that correspond to various sections of the teaching material. These scripts help students:

- Automate repetitive calculations and equation solving.
- Symbolically manipulate and analyze kinetic models efficiently.
- Visualize reaction dynamics and interpret chemical behavior using computational tools.

However, **all coding exercises in the "Hands-On Activities" sections** are **not** included in this repository, as students are encouraged to practice and develop their own solutions to build familiarity with coding.

## 📂 Repository Structure

Below is an overview of the folders included in this repository:

```
Chemical-Kinetics-Python
├── Chapter1
│   └── listing01.py
├── Chapter2
│   ├── listing01.py
│   ├── listing02.py
│   ├── listing03.py
│   ├── listing04.py
│   ├── listing05.py
│   ├── listing06.py
│   ├── listing07.py
│   ├── listing08.py
│   ├── listing09.py
│   └── listing10.py
├── Chapter3
│   ├── listing01.py
│   ├── listing02.py
│   ├── listing03.py
│   ├── listing04.py
│   ├── listing05.py
│   ├── listing06.py
│   ├── listing07.py
│   ├── listing08.py
│   └── listing09.py
├── Documents
│   └── book.pdf
├── LICENSE
└── README.md
```

- **Chapter1, Chapter2, Chapter3**: Contain Python scripts illustrating key concepts and examples for each chapter of the teaching material.
- **Documents**: Contains additional documentation or compiled references, e.g., a `book.pdf`.
- **LICENSE**: License information for this repository.
- **README.md**: The file you are currently reading.

## 🚀 Getting Started

### 1️⃣ Install Dependencies

Make sure you have Python and the required libraries installed:

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

Then, open or convert the `.py` files in the chapter folders into Jupyter Notebook cells for interactive exploration.

## 📌 How to Use the Code

Each script in the `ChapterX` folders corresponds to a section in the teaching material. You can:

- Run the provided scripts directly in a Python environment or in Jupyter Notebook.
- Modify them to see how different parameters affect reaction kinetics.
- Plot graphs or symbolic manipulations to deepen understanding of the reactions.

### Example: Simple Rate Equation

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

- **Practice the Hands-On Activities**: Attempt the coding exercises and derivations on your own, even though solutions are not provided here.
- **Experiment and Modify**: Change parameter values, add your own plots, or introduce new models to enhance your understanding.
- **Engage**: If you have questions, feel free to raise an issue on GitHub or discuss in class.

## 🐝 License

This repository is for educational purposes and is licensed under the terms specified in the [LICENSE](LICENSE) file.

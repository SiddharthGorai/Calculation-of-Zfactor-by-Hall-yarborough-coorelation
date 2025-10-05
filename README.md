# Z-Factor Calculator Using Hall-Yarborough Method

This project is a Python application that calculates the Z factor (compressibility factor) of a gas using the Hall-Yarborough method.
The Hall-Yarborough method is used in engineering to estimate the compressibility factor for gases, which is essential for various calculations in thermodynamics and fluid dynamics.
The application features a graphical user interface (GUI) built with Tkinter, allowing users to input necessary parameters and obtain the Z factor with ease.

### Features

- **Input Gas Composition:** Specify the composition of the gas mixture.
- **Reservoir Conditions:** Enter Pressure(Pa), Temperature(K), and Reduced density parameter.
- **Method Selection:** Choose among Kay's mixing rule, Standing's Correlation, and Sutton's Correlation for calculating critical pressure and temperature of the gas.
- **Hall-Yarborough Calculation:** Computes the Z-factor based on the input parameters and selected method.

## Requirements

```
numpy
matplotlib
tkinter 
```

## Installation

1. Clone or download this repository
2. Install required packages:
```bash
pip install numpy matplotlib
```

## Usage

### Running the Application

```bash
python main.py
```

### Input Parameters

1. **Pressure (Pa)**: Operating pressure in Pascals
2. **Temperature (K)**: Operating temperature in Kelvin
3. **Molecular Weight of Air (g/mol)**: Default value is 28.96 g/mol
4. **Reduced Density Parameter (Y₀)**: Initial guess for iteration (default: 0.01)
5. **Method for Ppr and Tpr Calculation**: Choose from three methods
6. **CSV Datafile**: Select composition data file

### CSV File Format

The composition CSV file should contain the following columns:

| Component | Molecular Weight (g/mol) | Mole Fraction | Critical Temp (°R) | Critical Pressure (psia) |
|-----------|-------------------------|---------------|-------------------|-------------------------|
| Methane   | 16.04                   | 0.85          | 343.0             | 667.8                   |
| Ethane    | 30.07                   | 0.10          | 549.8             | 707.8                   |
| ...       | ...                     | ...           | ...               | ...                     |

**Note**: The first row should be headers and will be ignored during calculations.

## Calculation Methods

### Method 1: Composition
Uses the composition data directly to calculate pseudo-critical properties:
- Ppc = Σ(yi × Pci)
- Tpc = Σ(yi × Tci)

### Method 2: Sutton's Correlation
Correlates pseudo-critical properties with specific gravity:
- Ppc = 756.8 - 131×SG - 3.6×SG²
- Tpc = 169.2 + 349.5×SG - 74×SG²

### Method 3: Standing's Correlation
Uses different correlations based on specific gravity:
- For SG < 0.75:
  - Ppc = 667 + 15×SG - 37.5×SG²
  - Tpc = 168 + 325×SG - 12.5×SG²
- For SG ≥ 0.75:
  - Ppc = 706 + 51.7×SG - 11.1×SG²
  - Tpc = 187 + 330×SG - 71.5×SG²

## Hall-Yarborough Method

The Hall-Yarborough method is an iterative approach to solve for the compressibility factor:

1. Calculate pseudo-reduced properties (Ppr, Tpr)
2. Set up Hall-Yarborough parameters (A, B, C, D)
3. Iterate using Newton-Raphson method until convergence
4. Calculate Z-factor from the converged reduced density

The iteration continues until |zi - z| < 10⁻⁵

## Output

The application provides:
- **Z-factor**: Gas compressibility factor
- **Gas Density**: In g/cm³

## References

- Hall, K.R., and Yarborough, L. "A New Equation of State for Z-Factor Calculations"


## Author

**Siddharth Gorai**  

- GitHub: [https://github.com/SiddharthGorai](https://github.com/SiddharthGorai)  
- LinkedIn: [https://www.linkedin.com/in/siddharth-gorai-ab01a7254/](https://www.linkedin.com/in/siddharth-gorai-ab01a7254/i)  
- Email: goraisiddharth@gmail.com

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
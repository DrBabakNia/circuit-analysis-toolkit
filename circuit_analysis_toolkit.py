import numpy as np
import cmath
import math
import re

class CircuitAnalysisToolkit:
    def __init__(self):
        self.student_name = ""
        
    def polar_to_rectangular(self, magnitude, angle_degrees):
        """Convert polar coordinates to rectangular form"""
        angle_radians = math.radians(angle_degrees)
        real = magnitude * math.cos(angle_radians)
        imag = magnitude * math.sin(angle_radians)
        return complex(real, imag)
    
    def rectangular_to_polar(self, complex_num):
        """Convert complex number to polar form"""
        magnitude = abs(complex_num)
        angle_rad = cmath.phase(complex_num)
        angle_deg = math.degrees(angle_rad)
        return magnitude, angle_deg
    
    def parse_complex_input(self, input_str):
        """
        Automatically parse complex number from various formats:
        - Rectangular: '10+20j', '10,20', '10 + 20j', '10'
        - Polar: '10<30', '10@30', '10 ∠ 30', '10 /_ 30'
        - Real: '10', '12.5'
        """
        input_str = input_str.strip().replace(' ', '')  # Remove spaces
        
        # If empty, return 0
        if not input_str:
            return complex(0)
        
        # Try to parse as real number first
        try:
            if 'j' not in input_str and '<' not in input_str and '@' not in input_str and '∠' not in input_str and '/_' not in input_str:
                return complex(float(input_str))
        except ValueError:
            pass
        
        # Pattern for rectangular format: a+bj or a,b
        rectangular_pattern = r'^([+-]?\d*\.?\d*)(?:([+-]\d*\.?\d*)j|,([+-]?\d*\.?\d*))?$'
        
        # Patterns for polar format
        polar_patterns = [
            r'^([+-]?\d*\.?\d*)[<@]([+-]?\d*\.?\d*)$',  # 10<30 or 10@30
            r'^([+-]?\d*\.?\d*)∠([+-]?\d*\.?\d*)$',     # 10∠30
            r'^([+-]?\d*\.?\d*)/_([+-]?\d*\.?\d*)$'     # 10/_30
        ]
        
        # Try rectangular format with 'j'
        match = re.match(rectangular_pattern, input_str)
        if match:
            real_part = match.group(1)
            imag_part_j = match.group(2)  # from a+bj format
            imag_part_comma = match.group(3)  # from a,b format
            
            try:
                real = float(real_part) if real_part not in ['', '+', '-'] else 0
                if real_part == '-':
                    real = -1.0
                
                if imag_part_j is not None:
                    imag = float(imag_part_j) if imag_part_j not in ['', '+', '-'] else 0
                    if imag_part_j == '-':
                        imag = -1.0
                    return complex(real, imag)
                elif imag_part_comma is not None:
                    imag = float(imag_part_comma)
                    return complex(real, imag)
                else:
                    return complex(real)  # Only real part
            except ValueError:
                pass
        
        # Try polar formats
        for pattern in polar_patterns:
            match = re.match(pattern, input_str)
            if match:
                try:
                    magnitude = float(match.group(1))
                    angle = float(match.group(2))
                    return self.polar_to_rectangular(magnitude, angle)
                except ValueError:
                    pass
        
        # If all parsing fails, try direct complex conversion
        try:
            # Handle cases like '2j', 'j', '-j'
            if input_str.endswith('j') and len(input_str) > 1:
                if input_str == 'j':
                    return complex(0, 1)
                elif input_str == '-j':
                    return complex(0, -1)
                else:
                    return complex(input_str)
            else:
                return complex(float(input_str))
        except ValueError:
            raise ValueError(f"Cannot parse input: {input_str}")
    
    def get_complex_input_flexible(self, prompt):
        """Get complex number input with automatic format detection"""
        print(f"\n{prompt}")
        print("Supported formats:")
        print("  Rectangular: '10+20j', '10,20', '10' (real only), '20j' (imaginary only)")
        print("  Polar: '10<30', '10@30', '10∠30', '10/_30'")
        
        while True:
            try:
                input_str = input("Enter value: ").strip()
                result = self.parse_complex_input(input_str)
                
                # Display what was recognized
                mag, angle = self.rectangular_to_polar(result)
                print(f"Recognized: {result.real:.4f} + {result.imag:.4f}j = {mag:.4f}∠{angle:.2f}°")
                
                return result
            except ValueError as e:
                print(f"Error: {e}")
                print("Please try again with supported format.")
            except Exception as e:
                print(f"Unexpected error: {e}")
                print("Please try again.")
    
    def display_complex_number(self, num, name="", unit=""):
        """Display complex number in both rectangular and polar formats"""
        magnitude, angle_deg = self.rectangular_to_polar(num)
        
        # Rectangular format
        if num.imag >= 0:
            rect_str = f"{num.real:.4f} + {num.imag:.4f}j"
        else:
            rect_str = f"{num.real:.4f} - {abs(num.imag):.4f}j"
        
        # Polar format
        polar_str = f"{magnitude:.4f} ∠ {angle_deg:.2f}°"
        
        if name:
            print(f"\n{name}:")
        print(f"  Rectangular: {rect_str} {unit}")
        print(f"  Polar: {polar_str} {unit}")
        
        return magnitude, angle_deg

    def display_equation_explanation(self):
        """Display explanation for equation solving"""
        print("\n" + "="*70)
        print("COMPLEX LINEAR EQUATIONS EXPLANATION")
        print("="*70)
        print("\nGeneral Form:")
        print("We solve systems of linear equations: A × x = b")
        print("\nWhere:")
        print("  A = Coefficient matrix (n×n complex numbers)")
        print("  x = Unknown vector (n×1 complex numbers)")
        print("  b = Result vector (n×1 complex numbers)")
        print("\nMatrix Form:")
        print("┌             ┐ ┌   ┐   ┌   ┐")
        print("│ A₁₁ A₁₂ A₁₃ │ │ x₁ │   │ b₁ │")
        print("│ A₂₁ A₂₂ A₂₃ │ × │ x₂ │ = │ b₂ │")
        print("│ A₃₁ A₃₂ A₃₃ │ │ x₃ │   │ b₃ │")
        print("└             ┘ └   ┘   └   ┘")
        print("\nApplications:")
        print("• Node voltage analysis in circuits")
        print("• Mesh current analysis")
        print("• AC circuit analysis with complex impedances")
        print("• Phasor domain circuit solutions")

    # =========================================================================
    # OPTION 1: SOLVE COMPLEX EQUATIONS
    # =========================================================================
    def solve_complex_equations(self):
        """Solve system of linear equations with complex coefficients"""
        self.display_equation_explanation()
        
        # Get matrix size
        while True:
            try:
                n = int(input("\nEnter number of equations (1-5): "))
                if 1 <= n <= 5:
                    break
                else:
                    print("Please enter a number between 1 and 5.")
            except ValueError:
                print("Please enter a valid integer.")
        
        # Initialize matrix and vector
        A = np.zeros((n, n), dtype=complex)
        b = np.zeros(n, dtype=complex)
        
        print(f"\nEnter coefficient matrix A ({n}x{n}):")
        for i in range(n):
            print(f"\nEquation {i+1} coefficients:")
            for j in range(n):
                prompt = f"  A[{i+1},{j+1}]:"
                A[i,j] = self.get_complex_input_flexible(prompt)
        
        print(f"\nEnter result vector b ({n}x1):")
        for i in range(n):
            prompt = f"  b[{i+1}]:"
            b[i] = self.get_complex_input_flexible(prompt)
        
        # Solve the system
        try:
            det = np.linalg.det(A)
            if abs(det) < 1e-10:
                print("Warning: Matrix is nearly singular! Solution may not be accurate.")
            
            x = np.linalg.solve(A, b)
            
            print("\n" + "="*60)
            print("SOLUTION RESULTS")
            print("="*60)
            
            for i, solution in enumerate(x):
                self.display_complex_number(solution, f"x[{i+1}]")
            
            print("\n" + "="*60)
            print("Analysis completed by: Dr. Babak Nia")
            print("Circuit Analysis Training Program")
            print("="*60)
                
        except np.linalg.LinAlgError:
            print("Error: Matrix is singular and cannot be solved.")

    def display_series_explanation(self):
        """Display explanation for series impedance"""
        print("\n" + "="*70)
        print("SERIES IMPEDANCE EXPLANATION")
        print("="*70)
        print("\nCircuit Diagram:")
        print("     Z₁      Z₂      Z₃           Zₙ")
        print("┌---███---███---███--- ... ---███---┐")
        print("│                                    │")
        print("V                                    │")
        print("│                                    │")
        print("└------------------------------------┘")
        print("\nFormula:")
        print("Z_eq = Z₁ + Z₂ + Z₃ + ... + Zₙ")
        print("\nProperties:")
        print("• Same current flows through all impedances")
        print("• Voltage divides across each impedance")
        print("• Total impedance is the sum of all impedances")
        print("• Phase angles add vectorially")

    # =========================================================================
    # OPTION 2: SERIES IMPEDANCE CALCULATION
    # =========================================================================
    def calculate_series_impedance(self):
        """Calculate equivalent impedance for series connection"""
        self.display_series_explanation()
        
        impedances = []
        
        while True:
            try:
                n = int(input("\nEnter number of impedances in series (1-10): "))
                if 1 <= n <= 10:
                    break
                else:
                    print("Please enter a number between 1 and 10.")
            except ValueError:
                print("Please enter a valid integer.")
        
        print(f"\nEnter {n} impedances:")
        for i in range(n):
            prompt = f"Impedance Z{i+1}"
            z = self.get_complex_input_flexible(prompt)
            impedances.append(z)
            self.display_complex_number(z, f"Z{i+1}", "Ω")
        
        # Calculate equivalent impedance
        z_eq = sum(impedances)
        
        print("\n" + "="*60)
        print("SERIES IMPEDANCE RESULTS")
        print("="*60)
        
        magnitude, angle_deg = self.display_complex_number(z_eq, "Equivalent Series Impedance", "Ω")
        
        # Additional analysis
        print(f"\nCircuit Analysis:")
        if z_eq.imag > 0:
            print("  Circuit Type: Inductive")
            print("  Power Factor: Lagging")
        elif z_eq.imag < 0:
            print("  Circuit Type: Capacitive") 
            print("  Power Factor: Leading")
        else:
            print("  Circuit Type: Purely Resistive")
            print("  Power Factor: Unity")
        
        power_factor = math.cos(math.radians(angle_deg))
        print(f"  Power Factor: {power_factor:.4f}")
        
        if z_eq.real != 0:
            Q = abs(z_eq.imag / z_eq.real)
            print(f"  Quality Factor (Q): {Q:.4f}")
        
        print("\n" + "="*60)
        print("Analysis completed by: Dr. Babak Nia")
        print("Circuit Analysis Training Program")
        print("="*60)

    def display_parallel_explanation(self):
        """Display explanation for parallel impedance"""
        print("\n" + "="*70)
        print("PARALLEL IMPEDANCE EXPLANATION")
        print("="*70)
        print("\nCircuit Diagram:")
        print("     ┌--- Z₁ ---┐")
        print("     ├--- Z₂ ---┤")
        print("  V  ├--- Z₃ ---┤")
        print("     ├   ...   ─┤")
        print("     └--- Zₙ ---┘")
        print("\nFormula:")
        print("1/Z_eq = 1/Z₁ + 1/Z₂ + 1/Z₃ + ... + 1/Zₙ")
        print("or using admittance: Y_eq = Y₁ + Y₂ + Y₃ + ... + Yₙ")
        print("\nProperties:")
        print("• Same voltage across all impedances")
        print("• Current divides among impedances")
        print("• Total admittance is the sum of all admittances")
        print("• Smaller impedance dominates the equivalent")

    # =========================================================================
    # OPTION 3: PARALLEL IMPEDANCE CALCULATION  
    # =========================================================================
    def calculate_parallel_impedance(self):
        """Calculate equivalent impedance for parallel connection"""
        self.display_parallel_explanation()
        
        impedances = []
        
        while True:
            try:
                n = int(input("\nEnter number of impedances in parallel (1-10): "))
                if 1 <= n <= 10:
                    break
                else:
                    print("Please enter a number between 1 and 10.")
            except ValueError:
                print("Please enter a valid integer.")
        
        print(f"\nEnter {n} impedances:")
        for i in range(n):
            prompt = f"Impedance Z{i+1}"
            z = self.get_complex_input_flexible(prompt)
            impedances.append(z)
            self.display_complex_number(z, f"Z{i+1}", "Ω")
        
        # Calculate equivalent impedance using admittances
        if n == 1:
            z_eq = impedances[0]
        else:
            y_total = sum(1/z for z in impedances)
            z_eq = 1 / y_total
        
        print("\n" + "="*60)
        print("PARALLEL IMPEDANCE RESULTS")
        print("="*60)
        
        magnitude, angle_deg = self.display_complex_number(z_eq, "Equivalent Parallel Impedance", "Ω")
        
        # Display admittance
        y_eq = 1 / z_eq
        self.display_complex_number(y_eq, "Equivalent Admittance", "S")
        
        # Additional analysis
        print(f"\nCircuit Analysis:")
        if z_eq.imag > 0:
            print("  Circuit Type: Inductive")
            print("  Power Factor: Lagging")
        elif z_eq.imag < 0:
            print("  Circuit Type: Capacitive")
            print("  Power Factor: Leading")
        else:
            print("  Circuit Type: Purely Resistive")
            print("  Power Factor: Unity")
        
        power_factor = math.cos(math.radians(angle_deg))
        print(f"  Power Factor: {power_factor:.4f}")
        
        print("\n" + "="*60)
        print("Analysis completed by: Dr. Babak Nia")
        print("Circuit Analysis Training Program")
        print("="*60)

    def display_delta_to_wye_explanation(self):
        """Display explanation for delta to wye conversion"""
        print("\n" + "="*70)
        print("DELTA TO WYE CONVERSION EXPLANATION")
        print("="*70)
        print("\nDelta Configuration (Δ):")
        print("        A")
        print("       / \\")
        print("   Z₁ /   \\ Z₃")
        print("     /     \\")
        print("    B---Z₂---C")
        print("\nWye Configuration (Y):")
        print("      A")
        print("      |")
        print("     Z_a")
        print("      |")
        print("B---Z_b---Z_c---C")
        print("      |")
        print("     Center")
        print("\nConversion Formulas:")
        print("Z_a = (Z₁·Z₃) / (Z₁ + Z₂ + Z₃)")
        print("Z_b = (Z₁·Z₂) / (Z₁ + Z₂ + Z₃)")
        print("Z_c = (Z₂·Z₃) / (Z₁ + Z₂ + Z₃)")
        print("\nApplications:")
        print("• Simplifying bridge circuits")
        print("• Three-phase system analysis")
        print("• Unbalanced load calculations")

    # =========================================================================
    # OPTION 4: DELTA TO WYE CONVERSION
    # =========================================================================
    def delta_to_wye(self):
        """Convert delta configuration to wye configuration"""
        self.display_delta_to_wye_explanation()
        
        print("\nEnter Delta impedances:")
        Z1 = self.get_complex_input_flexible("Delta impedance Z₁ (between nodes A-B)")
        Z2 = self.get_complex_input_flexible("Delta impedance Z₂ (between nodes B-C)")
        Z3 = self.get_complex_input_flexible("Delta impedance Z₃ (between nodes C-A)")
        
        denominator = Z1 + Z2 + Z3
        
        if abs(denominator) < 1e-12:
            print("Error: Sum of impedances is zero - conversion not possible")
            return
        
        # Calculate Wye impedances
        Za = (Z1 * Z3) / denominator  # Between node A and center
        Zb = (Z1 * Z2) / denominator  # Between node B and center  
        Zc = (Z2 * Z3) / denominator  # Between node C and center
        
        print("\n" + "="*60)
        print("DELTA TO WYE CONVERSION RESULTS")
        print("="*60)
        
        print("\nDelta Configuration:")
        self.display_complex_number(Z1, "Z₁ (A-B)", "Ω")
        self.display_complex_number(Z2, "Z₂ (B-C)", "Ω")
        self.display_complex_number(Z3, "Z₃ (C-A)", "Ω")
        
        print("\nWye Configuration:")
        self.display_complex_number(Za, "Z_a (A-center)", "Ω")
        self.display_complex_number(Zb, "Z_b (B-center)", "Ω")
        self.display_complex_number(Zc, "Z_c (C-center)", "Ω")
        
        print("\n" + "="*60)
        print("Analysis completed by: Dr. Babak Nia")
        print("Circuit Analysis Training Program")
        print("="*60)

    def display_wye_to_delta_explanation(self):
        """Display explanation for wye to delta conversion"""
        print("\n" + "="*70)
        print("WYE TO DELTA CONVERSION EXPLANATION")
        print("="*70)
        print("\nWye Configuration (Y):")
        print("      A")
        print("      |")
        print("     Z_a")
        print("      |")
        print("B---Z_b---Z_c---C")
        print("      |")
        print("     Center")
        print("\nDelta Configuration (Δ):")
        print("        A")
        print("       / \\")
        print("   Z₁ /   \\ Z₃")
        print("     /     \\")
        print("    B---Z₂---C")
        print("\nConversion Formulas:")
        print("Z₁ = (Z_a·Z_b + Z_b·Z_c + Z_c·Z_a) / Z_c")
        print("Z₂ = (Z_a·Z_b + Z_b·Z_c + Z_c·Z_a) / Z_a")
        print("Z₃ = (Z_a·Z_b + Z_b·Z_c + Z_c·Z_a) / Z_b")
        print("\nApplications:")
        print("• Circuit simplification")
        print("• Three-phase system transformations")
        print("• Network analysis")

    # =========================================================================
    # OPTION 5: WYE TO DELTA CONVERSION
    # =========================================================================
    def wye_to_delta(self):
        """Convert wye configuration to delta configuration"""
        self.display_wye_to_delta_explanation()
        
        print("\nEnter Wye impedances:")
        Za = self.get_complex_input_flexible("Wye impedance Z_a (node A to center)")
        Zb = self.get_complex_input_flexible("Wye impedance Z_b (node B to center)")
        Zc = self.get_complex_input_flexible("Wye impedance Z_c (node C to center)")
        
        numerator = Za*Zb + Zb*Zc + Zc*Za
        
        if abs(Zc) < 1e-12 or abs(Za) < 1e-12 or abs(Zb) < 1e-12:
            print("Error: Zero impedance - conversion not possible")
            return
        
        # Calculate Delta impedances
        Z1 = numerator / Zc  # Between nodes A-B
        Z2 = numerator / Za  # Between nodes B-C
        Z3 = numerator / Zb  # Between nodes C-A
        
        print("\n" + "="*60)
        print("WYE TO DELTA CONVERSION RESULTS")
        print("="*60)
        
        print("\nWye Configuration:")
        self.display_complex_number(Za, "Z_a (A-center)", "Ω")
        self.display_complex_number(Zb, "Z_b (B-center)", "Ω")
        self.display_complex_number(Zc, "Z_c (C-center)", "Ω")
        
        print("\nDelta Configuration:")
        self.display_complex_number(Z1, "Z₁ (A-B)", "Ω")
        self.display_complex_number(Z2, "Z₂ (B-C)", "Ω")
        self.display_complex_number(Z3, "Z₃ (C-A)", "Ω")
        
        print("\n" + "="*60)
        print("Analysis completed by: Dr. Babak Nia")
        print("Circuit Analysis Training Program")
        print("="*60)

def main():
    toolkit = CircuitAnalysisToolkit()
    
    print("\n" + "★" * 80)
    print("★" + " " * 78 + "★")
    print("★              CIRCUIT ANALYSIS TOOLKIT - COMPREHENSIVE SUITE                ★")
    print("★" + " " * 78 + "★")
    print("★                  Instructor: Dr. Babak Nia                                ★")
    print("★                  Course: Circuit Analysis Training                        ★")
    print("★" + " " * 78 + "★")
    print("★" * 80)
    
    toolkit.student_name = input("\nEnter your name: ")
    
    while True:
        print(f"\nWelcome {toolkit.student_name}!")
        print("\n" + "="*70)
        print("MAIN MENU - CIRCUIT ANALYSIS TOOLKIT")
        print("="*70)
        print("1. Solve Complex Linear Equations (up to 5×5)")
        print("2. Calculate Series Impedance (up to 10 impedances)")
        print("3. Calculate Parallel Impedance (up to 10 impedances)")
        print("4. Delta to Wye Conversion (Δ → Y)")
        print("5. Wye to Delta Conversion (Y → Δ)")
        print("6. Input Format Examples")
        print("7. Exit")
        
        choice = input("\nEnter your choice (1-7): ").strip()
        
        if choice == '1':
            toolkit.solve_complex_equations()
        elif choice == '2':
            toolkit.calculate_series_impedance()
        elif choice == '3':
            toolkit.calculate_parallel_impedance()
        elif choice == '4':
            toolkit.delta_to_wye()
        elif choice == '5':
            toolkit.wye_to_delta()
        elif choice == '6':
            show_input_examples(toolkit)
        elif choice == '7':
            print(f"\nThank you for using the Circuit Analysis Toolkit, {toolkit.student_name}!")
            print("Good luck with your circuit analysis studies!")
            break
        else:
            print("Invalid choice! Please enter a number between 1 and 7.")
        
        input("\nPress Enter to continue...")

def show_input_examples(toolkit):
    """Show examples of supported input formats"""
    print("\n" + "="*70)
    print("INPUT FORMAT EXAMPLES")
    print("="*70)
    
    examples = [
        ("10+20j", "Rectangular: 10 + 20j"),
        ("10,20", "Rectangular: 10 + 20j (comma separated)"),
        ("10", "Real number: 10 + 0j"),
        ("20j", "Imaginary only: 0 + 20j"),
        ("-15+30j", "Negative real: -15 + 30j"),
        ("10-5j", "Negative imaginary: 10 - 5j"),
        ("10<30", "Polar: 10∠30° (magnitude < angle)"),
        ("10@30", "Polar: 10∠30° (magnitude @ angle)"),
        ("10∠30", "Polar: 10∠30° (magnitude ∠ angle)"),
        ("10/_30", "Polar: 10∠30° (magnitude /_ angle)"),
        ("5.5+2.3j", "Decimal: 5.5 + 2.3j"),
        ("100", "Large real: 100 + 0j"),
    ]
    
    print("\nTry these examples in any input field:\n")
    for input_example, description in examples:
        try:
            result = toolkit.parse_complex_input(input_example)
            mag, angle = toolkit.rectangular_to_polar(result)
            print(f"  Input: '{input_example:12}' → {description:45} → Result: {result.real:6.2f} + {result.imag:6.2f}j = {mag:6.2f}∠{angle:5.1f}°")
        except Exception as e:
            print(f"  Input: '{input_example:12}' → Error: {e}")

if __name__ == "__main__":
    main()
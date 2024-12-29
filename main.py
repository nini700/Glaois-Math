class GaloisField:
    """
    A class implementing Galois Field arithmetic operations over GF(2^m).
    This implementation uses polynomial representation with coefficients in GF(2).
    """
    def __init__(self, primitive_poly, m):
        """
        Initialize Galois Field with primitive polynomial and field degree.
        
        Args:
            primitive_poly (int): The primitive polynomial in binary representation
            m (int): The degree of the field (GF(2^m))
        """
        self.m = m
        self.primitive_poly = primitive_poly
        self.field_size = 1 << m
        self.mask = self.field_size - 1
        
        self.exp_table = [0] * (2 * self.field_size)
        self.log_table = [0] * self.field_size
        
        x = 1
        for i in range(self.field_size - 1):
            self.exp_table[i] = x
            self.exp_table[i + self.field_size - 1] = x
            self.log_table[x] = i
            x = self._multiply_without_tables(x, 2)
    
    def _multiply_without_tables(self, a, b):
        """Helper function for multiplication without using lookup tables."""
        result = 0
        while b:
            if b & 1:
                result ^= a
            a <<= 1
            if a & self.field_size:
                a ^= self.primitive_poly
            b >>= 1
        return result
    
    def add(self, a, b):
        """Addition in GF(2^m) is just XOR."""
        return a ^ b
    
    def subtract(self, a, b):
        """Subtraction in GF(2^m) is the same as addition."""
        return self.add(a, b)
    
    def multiply(self, a, b):
        """Multiplication in GF(2^m) using lookup tables."""
        if a == 0 or b == 0:
            return 0
        return self.exp_table[(self.log_table[a] + self.log_table[b]) % (self.field_size - 1)]
    
    def divide(self, a, b):
        """Division in GF(2^m) using lookup tables."""
        if b == 0:
            raise ValueError("Division by zero")
        if a == 0:
            return 0
        return self.exp_table[(self.log_table[a] - self.log_table[b]) % (self.field_size - 1)]
    
    def inverse(self, a):
        """Find multiplicative inverse in GF(2^m)."""
        if a == 0:
            raise ValueError("Zero has no multiplicative inverse")
        return self.exp_table[(-self.log_table[a]) % (self.field_size - 1)]

# Example usage
def example_usage():
    gf = GaloisField(primitive_poly=0b10011, m=4)
    a, b = 0b1011, 0b1101  # 11 and 13 in decimal
    
    print(f"Field: GF(2^4)")
    print(f"a = {bin(a)[2:].zfill(4)}, b = {bin(b)[2:].zfill(4)}")
    print(f"a + b = {bin(gf.add(a, b))[2:].zfill(4)}")
    print(f"a * b = {bin(gf.multiply(a, b))[2:].zfill(4)}")
    print(f"a / b = {bin(gf.divide(a, b))[2:].zfill(4)}")
    print(f"inverse(a) = {bin(gf.inverse(a))[2:].zfill(4)}")

if __name__ == "__main__":
    example_usage()
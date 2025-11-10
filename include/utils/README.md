# File Reading Utilities

This directory contains utility functions for reading data from files.

## Files

### `readinteger.h`
- **Function**: `readinteger(FILE *file)`
- **Description**: Reads a positive integer from a file
- **Returns**: The integer read, or `INT32_MAX` if EOF is reached
- **Usage**: For reading identifiers, counters, etc.

### `readreal.h`
- **Function**: `readreal(FILE *file)`
- **Description**: Reads a positive real number from a file
- **Support**: Decimal notation (e.g., 123.456)
- **Returns**: The floating point number read, or `FLT_MAX` if EOF is reached
- **Usage**: For reading frequencies, probabilities, etc.

### `readnegativereal.h`
- **Function**: `readnegativereal(FILE *file)`
- **Description**: Reads a real number (positive or negative) from a file
- **Support**: 
  - Negative numbers (with '-' prefix)
  - Decimal notation (e.g., -123.456)
  - Scientific notation (e.g., 1.23e-4)
- **Returns**: The double precision number read, or `FLT_MAX` if EOF is reached
- **Usage**: For reading signed values, correlations, etc.

### `FileIOUtils.h`
- **Description**: Convenience header including all utilities
- **Usage**: `#include "utils/FileIOUtils.h"` to include all functions

## Usage Example

```cpp
#include "utils/FileIOUtils.h"
#include <stdio.h>

FILE* file = fopen("data.txt", "r");
if (file) {
    int id = readinteger(file);
    float frequency = readreal(file);
    double correlation = readnegativereal(file);
    fclose(file);
}
```

## Notes

- All functions automatically skip non-numeric characters
- Functions return maximum values (INT32_MAX, FLT_MAX) on EOF
- Functions are inline and defined in headers for optimization

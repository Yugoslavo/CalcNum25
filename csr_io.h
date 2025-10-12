#ifndef CSR_IO_H
#define CSR_IO_H

#ifdef __cplusplus
extern "C" {
#endif

// Writes CSR arrays to files: <prefix>_ia.bin, <prefix>_ja.bin, <prefix>_a.bin
// Return 0 on success, non-zero on error.
int write_csr_arrays(const char *prefix, int n,
                     const int *ia, const int *ja, const double *a);

#ifdef __cplusplus
}
#endif
#endif // CSR_IO_H
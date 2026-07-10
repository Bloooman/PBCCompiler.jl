// Minimal subset of the OpenQASM 2.0 standard header (qelib1.inc)
// containing the gates needed by the test fixtures. Definitions follow
// the Quantum Experience standard library.

// Toffoli
gate ccx a,b,c
{
  h c;
  cx b,c; tdg c;
  cx a,c; t c;
  cx b,c; tdg c;
  cx a,c; t b; t c; h c;
  cx a,b; t a; tdg b;
  cx a,b;
}

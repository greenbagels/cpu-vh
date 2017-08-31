//return a 4 dimensional linear interpolation inside the hypercube, given the values
//on the corners (a0000 through a1111) and edge lengths x0 through x3
double linearInterp4D(double x0, double x1, double x2, double x3,
                      double a0000, double a1000, double a0100, double a0010, double a0001,
                      double a1100, double a1010, double a1001,
                      double a0110, double a0101, double a0011,
                      double a1110, double a1101, double a0111, double a1011, double a1111)
{
  double result = 0;
  result = ((1-x0) * (1-x1) * (1-x2) * (1-x3) * a0000)
            + ((x0) * (1-x1) * (1-x2) * (1-x3) * a1000)
            + ((1-x0) * (x1) * (1-x2) * (1-x3) * a0100)
            + ((1-x0) * (1-x1) * (x2) * (1-x3) * a0010)
            + ((1-x0) * (1-x1) * (1-x2) * (x3) * a0001)
            + ((x0) * (x1) * (1-x2) * (1-x3) * a1100)
            + ((x0) * (1-x1) * (x2) * (1-x3) * a1010)
            + ((x0) * (1-x1) * (1-x2) * (x3) * a1001)
            + ((1-x0) * (x1) * (x2) * (1-x3) * a0110)
            + ((1-x0) * (x1) * (1-x2) * (x3) * a0101)
            + ((1-x0) * (1-x1) * (x2) * (x3) * a0011)
            + ((x0) * (x1) * (x2) * (1-x3) * a1110)
            + ((x0) * (x1) * (1-x2) * (x3) * a1101)
            + ((x0) * (1-x1) * (x2) * (x3) * a1011)
            + ((1-x0) * (x1) * (x2) * (x3) * a0111)
            + ((x0) * (x1) * (x2) * (x3) * a1111);

  return result;
}

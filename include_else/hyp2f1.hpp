#ifndef HYP2F1_INCLUDED
#define HYP2F1_INCLUDED




using namespace std;

#define SIGN(a) (((a) < 0) ? (-1) : (1))

double hyp2f1a(double beta, double z);
double hyp2f1b(double beta, double z);
double hyp2f1c(double gamma, double z);
double hyp2f1d(double gamma, double z);

#endif // HYP2F1_INCLUDED
#ifndef DERIVATIVES_H
#define DERIVATIVES_H

struct params;
extern int g_NSYS;
extern struct params g_consts;
void set_parameters(struct params *);
void derivative(const double *, const double [], double []);
void derivative_internal(const double *, const double [], double [], const struct params);
int events(const double *, const double []);

#endif // DERIVATIVES_H

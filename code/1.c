#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_rng.h>
/////////////////////////////////////
int main(void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for (i = 0; i < n; i++) 
    {
      double u = gsl_rng_uniform (r);
      printf ("%.5f\n", u);
    }

  printf ("generator type: %s\n", gsl_rng_name (r));
  printf ("seed = %lu\n", gsl_rng_default_seed);
  printf ("first value = %lu\n", gsl_rng_get (r));


  gsl_rng_free (r);

  return 0;
}

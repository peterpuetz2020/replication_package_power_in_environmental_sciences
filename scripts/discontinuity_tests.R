## Backward-compatible entry point. The maintained implementation is caliper_tests.R.
source(here::here("scripts", "caliper_tests.R"))
if (sys.nframe() == 0L) run_discontinuity_tests()

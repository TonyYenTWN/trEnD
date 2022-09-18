// Source file for power flow analysis
#include "power_network.h"
#include "src/power_market/power_market.h"

// HELM method
// P-Q Bus
// Set {V}(s) = \sum {a}_n * s^n,
// {\hat V}(s) = \sum {b}_n * s^n,
// {1 / V}(s) = \sum {c}_n * s^n,
// {1 / \hat V}(s) = \sum {d}_n * s^n,

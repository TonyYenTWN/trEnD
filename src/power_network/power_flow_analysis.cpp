// Source file for power flow analysis
#include "power_network.h"
#include "src/power_market/power_market.h"

// HELM method
// P-Q buses
// Set {V}(s) = \sum {a}_n * s^n
// {\hat V}(s) = \sum {b}_n * s^n
// {1. / V}(s) = \sum {c}_n * s^n
// {1. / \hat V}(s) = \sum {d}_n * s^n
// Equations:
// [Y] {V}(s) = s * {Conj(S)} . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = s * {S} . {1. / V}(s)
//
// P-U buses
// Set {a}, {b}, {c}, {d} the same as P-Q bus
// {Q}(s) = \sum {e}_n * s^n
// Equations:
// [Y] {V}(s) = s * ({P} - j {Q}(s)) . {1. / \hat V}(s)
// [Conj(Y)] {\hat V}(s) = s * ({P} + j {Q}(s)) . {1. / V}(s)
// {V}(s) * {\hat V}(s) = {U_0^2} + s * ({U^2} - {U_0^2})
//
// Conservation law for currents
// {1}^T [Y] {V}(s) = 0.


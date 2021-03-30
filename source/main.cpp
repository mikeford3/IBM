#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

// https://www.rand.org/content/dam/rand/pubs/papers/2008/P4466.pdf
// https://apps.dtic.mil/sti/pdfs/ADA399211.pdf

const auto R  = 8.31441; // J / (K * mol), has constant
const auto ra = 0.69;    // relative area of propellant cylinders to the total area of the gun/shell

struct Propellant {
    double mass;          // kg
    double mass_0;        // kg
    double density;       // kg / m^3
    double area;          // total area of the cylinders in the propellant (less than the area of the barrel)
    double a;             // burn rate coefficient, m/s
    double n;             // burn rate exponent
    double length;        // metres
    double energyDensity; // J / kg
};

struct Gas {
    double covolume; // m^3 / kg - the volume of gas created from 1 kg of solid;
    double n;        // # number of mols
    double molD;     // # / kg - number of mols per kg
    double P;        // Pa, pressure
    double T;        // Kelvin, temperature
    double u;        // J, energy
    double u0;       // J, initial energy

    double specHeat; // J / (mol K)
};

struct Shell {
    double mass; // kg
    double area; // of base, m^2
    // TODO - change the frictions to a friction coefficient?
    double frictionS; // static friction force, N
    double frictionD; // dynamic friction, N / (m/s)
    double speed;     // m / s
};

struct Gun {
    double volumeC; // m^3, free chamber volume (i.e. space left by charge burnup & space between cylinders in charge)
                    // plus change in volume from shell moving
    double length;  // m, barell length
    double d;       // m, position of base of shell
    double area;    // same as shell, m^2
};

struct Force {
    double fg; // force applied by the gas (i.e. excluding friction etc and used to calculate work done by gas)
    double fs; // net force applied to shell after friction removed etc
};

/** Burn propellant, add matter and energy to gas whilst reducing mass of propellant,
 *  calculate updated temperature, substance and pressure of gas
 *
 *  'p' and 'g' are in/out args
 *  'dt' is in arg, time step
 */
void burn(Propellant& p, Gas& g, double dt) noexcept {
    auto mb = p.a * std::pow(g.P, p.n) * p.area * dt * p.mass / p.mass_0; // mass burned
    if (mb <= p.mass) {
        p.mass -= mb;
    } else {
        mb     = p.mass;
        p.mass = 0;
    }
    auto du = mb * p.energyDensity; // additional energy added
    auto dn = mb * g.molD;          // additional substance added

    g.P = g.P * ((g.u + du) / g.u) * (g.n + dn) / g.n; // new pressure, increase with energy & n

    g.u += du;
    g.n += dn;
}

/** Calculate force applied to shell based on gas pressure
 *  's' is the shell, in arg
 *  'g' is the gas, in arg
 *  Returns the forces N, cal be -ve
 *
 */
Force calculate_force(const Shell& s, const Gas& g) noexcept {
    auto f = Force{};
    f.fg   = g.P * s.area;
    f.fs   = (g.P - 101e3) * s.area; // net force of gas on base
    if (s.speed == 0) {
        // shell not moving, check against static friction, return 0 if not overcome
        // otherwise return the surplus force to get the shell moving
        f.fs = std::max(f.fs - s.frictionS, 0.);
    } else {
        // shell moving, remove dynamic friction force
        f.fs -= s.frictionD * s.speed;
    }
    return f;
}

/** Update gas due to shell movement
 *  's' is the shell, in arg
 *  'ga' is the gas, in/out arg
 *  'gu' is the gun, in/out arg, volume updated
 *  'f'  force applied to the shell/by the gas
 *  'd_0' is the old shell distance
 *  'dd' is the distance the base of the shell has moved
 *  Remove energy in gas due to work done on shell
 *  Increase the volume of the gas due to shell movement down barrell
 *  Reduce pressure of gas due to volume increase & loss of energy
 *
 */
void expand_gas(const Shell& s, Gas& ga, Gun& gu, Force f, double dd) noexcept {
    // account for movement of shell on volume and pressure
    auto dv = dd * gu.area; // extra volume due to shell movement

    // auto pv_0 = ga.P * gu.volumeC;                     // old pressure * volume
    ga.P = ga.P * gu.volumeC / (gu.volumeC + dv); // update pressure, accounting for change in volume only
    gu.volumeC += dv;                             // new volume

    // account for work done by gas, update pressure and temperature
    auto du = f.fg * dd;        // work done by gas to shell, always positive (shell not allowed to move backwards)
    ga.P *= (ga.u - du) / ga.u; // update pressure of gas in in same way as temperature
    ga.u -= du;
    ga.P = std::max(ga.P, 20.e3); // arbitrary
    ga.u = std::max(ga.u, ga.u0);
}

int main() {
    // 15in shell
    auto s = Shell{.mass      = 879,
                   .area      = std::numbers::pi * 0.19 * 0.19,
                   .frictionS = 87900, // shell moves once approx 10x its weight is taken up (wild guess)
                   .frictionD = 10000, //
                   .speed     = 0.};

    auto p = Propellant{
        .mass          = 196,
        .mass_0        = 196,
        .density       = 1.5e3,
        .area          = s.area * ra,
        .a             = 8000, // made up!
        .n             = 0.15, // made up!
        .length        = 196 / (1.5e3 * s.area * ra),
        .energyDensity = 3349'440 // 800 cal/gram in J / kg (then div/10 to give more sensitive answer)
    };

    auto gu = Gun{.volumeC = 0.502 * (1. - ra), // m^3 30,650 in^3
                  .length  = .38 * 42,
                  .d       = 0.,
                  .area    = s.area};

    auto ga = Gas{
        .covolume = 9.5e-4,                      // 9.5 cm^3/gram in m^3/kg, nitrocellulose
                                                 // (not used anywhere, probably due to universal gas law)
        .n        = gu.volumeC * 1.292 * 0.04e3, // vol of chamber [m^3] * density of air [kg/m^3] / [mol / kg]
        .molD     = 0.04e3,                      // 0.04 mol / gram
        .P        = 101e3,
        .T        = 293.15,
        .u        = 293.15 * 36.94 * gu.volumeC * 1.292, //* 0.04e3,        // to T * specHeat * n
        .u0       = 293.15 * 36.94 * gu.volumeC * 1.292, //* 0.04e3,        // to T * specHeat * n
        .specHeat = 36.94                                // J / (mol K), CO2 (don't know for product gases)
    };

    printf("Starting analysis of 15in shell\n");
    printf("%10s %10s %10s %10s %10s %10s %10s %10s\n", "step", "time", "powder", "pressure", "temperature", "volume",
           "distance", "speed");
    printf("%10s %10s %10s %10s %10s %10s %10s %10s\n", "#", "s", "kg", "kPa", "degC", "m^3", "m", "m/s");

    const auto dt = 0.000001; // time step, 1 ms

    for (auto i = 0; i < 10000000; ++i) {
        // burn propellant, create gas and add energy to gas
        burn(p, ga, dt);

        // gas applies force to shell
        auto f = calculate_force(s, ga);

        // shell accelerate
        auto a = f.fs / s.mass;
        s.speed += a * dt;

        // shell moves
        auto dd = s.speed * dt;
        gu.d += dd; // implicit euler
        gu.d = std::min(gu.d, gu.length);

        // update gas,
        // reduce energy in gas due to work on shell (should balance with KE increase of shell)

        // volume of gas increases to take up space of shell, reduce pressure
        expand_gas(s, ga, gu, f, dd);
        if (i % 1000 == 0)
            printf("%10i %10.6f %10.2f %10.0f %10.1f %10.3f %10.3f %10.2f\n", i, i * dt, p.mass, ga.P / 1000,
                   ga.T - 273.15, gu.volumeC, gu.d, s.speed);
        // shell left gun, exit
        if (gu.d >= gu.length) {
            if (i % 1000 != 0)
                printf("%10i %10.6f %10.2f %10.0f %10.1f %10.3f %10.3f %10.2f\n", i, i * dt, p.mass, ga.P / 1000,
                       ga.T - 273.15, gu.volumeC, gu.d, s.speed);
            break;
        }
    }
}
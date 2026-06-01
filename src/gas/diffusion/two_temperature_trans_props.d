/**
 * two_temperature_trans_props.d
 *
 * Author: Rowan G.
 * Date: 2021-03-16
 * History: 2021-03-16 -- code ported from two_temperature air model.
 *
 * References
 * ----------
 * Gnoffo, Gupta and Shin (1989)
 * Conservation equations and physical models for hypersonic air flows
 * in thermal and chemical nonequilibrium.
 * NASA TP-2867, February 1989, NASA Langley Research Center
 *
 * Gupta, Yos and Thompson (1989)
 * A review of reaction rate and thermodynamic and transport properties
 * for the 11-species air model for chemical and thermal nonequilibrium
 * calculations to 30,000 K.
 * NASA TM-101528, February 1989, NASA Langley Research Center
 **/

module gas.diffusion.two_temperature_trans_props;

import std.math;
import std.string;
import std.conv;
import std.algorithm;

import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas;
import gas.diffusion.transport_properties_model;
import gas.physical_constants;

enum CUSTOM_LAMP_N_SPECIES = 11;
enum CUSTOM_LAMP_HEAVY_SPECIES = 10;
enum CUSTOM_LAMP_THERMAL_SYSTEM = 20;
enum CUSTOM_LAMP_AVOGADRO_NUMBER = 6.02214076e23;
enum double CUSTOM_LAMP_PI = 3.141592653589793238462643383279502884;
enum double CUSTOM_LAMP_CS = 1.0e-20;
enum double CUSTOM_LAMP_KB = 1.380649e-23;
enum double CUSTOM_LAMP_EC = 1.602176634e-19;
enum double CUSTOM_LAMP_EPS0 = 8.8541878188e-12;
enum double CUSTOM_LAMP_EULER_GAMMA = 0.5772156649015329;
enum double CUSTOM_LAMP_H = 6.62607015e-34;
enum double CUSTOM_LAMP_C0 = 299792458.0;
enum double CUSTOM_LAMP_EV2J = CUSTOM_LAMP_EC;
enum double CUSTOM_LAMP_AMU = 1.66053906892e-27;
enum double CUSTOM_LAMP_ME = 9.1093837139e-31;
enum double CUSTOM_LAMP_HVIC_DENOM_FLOOR = 1.0e-200;
enum double CUSTOM_LAMP_PIVOT_FLOOR = 1.0e-300;
enum double CUSTOM_LAMP_INTERNAL_DENOM_FLOOR = 1.0e-300;

@nogc number customLAMPNumberFromDouble(double value)
{
    version(complex_numbers) {
        return complex(value);
    } else {
        return value;
    }
}

@nogc double customLAMPReal(in number value)
{
    version(complex_numbers) {
        return value.re;
    } else {
        return value;
    }
}

class TwoTemperatureTransProps : TransportPropertiesModel {
public:

    this(lua_State *L, string[] speciesNames)
    {
        mNSpecies = to!int(speciesNames.length);
        double[] molMasses;
        mMolMasses.length = mNSpecies;
        mParticleMass.length = mNSpecies;
        mMolef.length = mNSpecies;
        mCharge.length = mNSpecies;
        allocateArraysForCollisionIntegrals();

        foreach (isp, spName; speciesNames) {
            if (spName == "e-") mElectronIdx = to!int(isp);
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);
            double m = getDouble(L, -1, "M");
            mMolMasses[isp] = m;
            mParticleMass[isp] = 1000.0*m/Avogadro_number; // kg -> g
            mCharge[isp] = getInt(L, -1, "charge");
            string type = getString(L, -1, "type");
            if (type == "molecule") {
                mMolecularSpecies ~= to!int(isp);
            }
            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        // Fill in support data for collision integrals
        lua_getglobal(L, "db");
        lua_getfield(L, -1, "CIs");
        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                string key = speciesNames[isp] ~ ":" ~ speciesNames[jsp];
                lua_getfield(L, -1, key.toStringz);
                if (lua_isnil(L, -1)) {
                    lua_pop(L, 1);
                    // Try reverse order, eg. N2:O2 --> O2:N2
                    key = speciesNames[jsp] ~ ":" ~ speciesNames[isp];
                    lua_getfield(L, -1, key.toStringz);
                    if (lua_isnil(L, -1)) {
                        // There's still a problem, can't find entry for CI data.
                        string msg = format("Collision integral data for '%s' is missing.\n", key);
                        throw new Error(msg);
                    }
                }
                lua_getfield(L, -1, "pi_Omega_11");
                mA_11[isp][jsp] = getDouble(L, -1, "A");
                mB_11[isp][jsp] = getDouble(L, -1, "B");
                mC_11[isp][jsp] = getDouble(L, -1, "C");
                mD_11[isp][jsp] = getDouble(L, -1, "D");
                lua_pop(L, 1); // pop: pi_Omega_11
                lua_getfield(L, -1, "pi_Omega_22");
                mA_22[isp][jsp] = getDouble(L, -1, "A");
                mB_22[isp][jsp] = getDouble(L, -1, "B");
                mC_22[isp][jsp] = getDouble(L, -1, "C");
                mD_22[isp][jsp] = getDouble(L, -1, "D");
                lua_pop(L, 1); // pop: pi_Omega_22
                lua_pop(L, 1); // pop: collision pair.
            }
        }
        lua_pop(L, 1); // pop: cis
        lua_pop(L, 1); // pop: db

        // Compute alphas
        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. mNSpecies) {
                double M_isp = mMolMasses[isp];
                double M_jsp = mMolMasses[jsp];
                mMu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                mMu[isp][jsp] *= 1000.0; // convert kg/mole --> g/mole
                double M_ratio = M_isp/M_jsp;
                double numer = (1.0 - M_ratio)*(0.45 - 2.54*M_ratio);
                double denom = (1.0 + M_ratio)^^2;
                mAlpha[isp][jsp] = 1.0 + numer/denom;
                mIsCoulombCollision[isp][jsp] = (mCharge[isp]!=0)&&(mCharge[jsp]!=0);
            }
        }
    }

    void allocateArraysForCollisionIntegrals()
    {
        mMu.length = mNSpecies;
        mAlpha.length = mNSpecies;
        mA_11.length = mNSpecies;
        mB_11.length = mNSpecies;
        mC_11.length = mNSpecies;
        mD_11.length = mNSpecies;
        mDelta_11.length = mNSpecies;
        mA_22.length = mNSpecies;
        mB_22.length = mNSpecies;
        mC_22.length = mNSpecies;
        mD_22.length = mNSpecies;
        mDelta_22.length = mNSpecies;
        mIsCoulombCollision.length = mNSpecies;
        foreach (isp; 0 .. mNSpecies) {
            mMu[isp].length = mNSpecies;
            mAlpha[isp].length = mNSpecies;
            mA_11[isp].length = isp+1;
            mB_11[isp].length = isp+1;
            mC_11[isp].length = isp+1;
            mD_11[isp].length = isp+1;
            mDelta_11[isp].length = mNSpecies;
            mA_22[isp].length = isp+1;
            mB_22[isp].length = isp+1;
            mC_22[isp].length = isp+1;
            mD_22[isp].length = isp+1;
            mDelta_22[isp].length = mNSpecies;
            mIsCoulombCollision[isp].length = mNSpecies;
        }
    }

    @nogc
    override void updateTransProps(ref GasState gs, GasModel gm)
    {
        massf2molef(gs.massf, mMolMasses, mMolef);
        // Computation of transport coefficients via collision integrals.
        // Equations follow those in Gupta et al. (1990)
        computeDelta22(gs);
        computeDelta11(gs);

        // Compute mixture viscosity.
        number sumA = 0.0;
        number denom;
        foreach (isp; 0 .. mNSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_22[isp][jsp];
            }
            if (isp == mElectronIdx) continue;
            sumA += mParticleMass[isp]*mMolef[isp]/denom;
        }
        // Add term if electron present.
        if (mElectronIdx != -1) {
            // An additional term is required in the mixture viscosity.
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_22[mElectronIdx][jsp];
            }
            sumA += mParticleMass[mElectronIdx]*mMolef[mElectronIdx]/denom;

        }
        gs.mu = sumA * (1.0e-3/1.0e-2); // convert g/(cm.s) -> kg/(m.s)

        // Compute component thermal conductivities
        // k in transrotational = k_tr + k_rot
        // k in vibroelectronic = k_ve + k_E
        // 1. k_tr
        sumA = 0.0;
        foreach (isp; 0 .. mNSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                if (jsp != mElectronIdx) {
                    denom += mAlpha[isp][jsp]*mMolef[jsp]*mDelta_22[isp][jsp];
                }
                else {
                    denom += 3.54*mAlpha[isp][jsp]*mMolef[jsp]*mDelta_22[isp][jsp];
                }
            }
            if (isp == mElectronIdx) continue;
            sumA += mMolef[isp]/denom;
        }
        double kB_erg = 1.38066e-16; // erg/K
        number k_tr = 2.3901e-8*(15./4.)*kB_erg*sumA;
        k_tr *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)

        // 2. k_rot and k_vib at the same time
        // Assuming rotation is fully excited, eq (75) in Gnoffo
        // and vibration is partially excited eq (52a) in Gupta, Yos, Thomson
        number k_rot = 0.0;
        number k_vib = 0.0;

        // For k_vib, Cp(T_vib) needs to be evaluated
        foreach (isp; mMolecularSpecies) {
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += mMolef[jsp]*mDelta_11[isp][jsp];
            }
            k_rot += mMolef[isp]/denom;

            // figure out how much vibrational excitation there is
            // The equation in the paper doesn't call for the min,
            // but it has been included to get the same result as
            // when the flow is fully excited
            number Cp = mMolMasses[isp]/R_universal*gm.Cp(gs, isp);
            k_vib += fmax((fmin(Cp.re, 9.0/2.0) - 7.0/2.0), 0.0) * mMolef[isp]/denom;
        }
        k_rot *= 2.3901e-8*kB_erg;
        k_rot *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        k_vib *= 2.3901e-8*kB_erg;
        k_vib *= (4.184/1.0e-2); // cal/cm.s.K) --> J/(m.s.K)

        // Eq (76) in Gnoffo.
        gs.k = k_tr + k_rot;

        // 4. k_e
        number k_E = 0.0;
        if (mElectronIdx != -1) {
            // electron present.
            denom = 0.0;
            foreach (jsp; 0 .. mNSpecies) {
                denom += 1.45*mMolef[jsp]*mDelta_22[mElectronIdx][jsp];
            }
            k_E = mMolef[mElectronIdx]/denom;
            k_E *= 2.3901e-8*(15./4.)*kB_erg;
            k_E *= (4.184/1.0e-2); // cal/(cm.s.K) --> J/(m.s.K)
        }
        gs.k_modes[0] = k_vib + k_E;
    }

    @nogc void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D)
    {
    /*
        Compute the binary diffusion coefficients using the Omega_11 collisions
        integrals and equations 42a from Gupta et al. 1990.

        Notes:
         - Prior to December 2021, this routine was returning diffusion coefficients
           that were 100 times too small. This was due to a unit conversion error
           related to Delta_11, which is in cm-s, not m-s. (NNG)

         - TODO: Think about "p" in this equation. Should it be total pressure, or
           bath pressure of binary interactors only? (RJG)
    */
        massf2molef(gs.massf, mMolMasses, mMolef);
        computeDelta11(gs);
        immutable double kB_erg = 1.38066e-16; // erg/K
        number p_cgs = gs.p*10.0; // CGS pressure in Baryes (?!)

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                number T = (isp == mElectronIdx || jsp == mElectronIdx) ? gs.T_modes[0] : gs.T;
                number Dij = (kB_erg * T)/(p_cgs * mDelta_11[isp][jsp]); // cm^2/s
                Dij /= (100*100); // cm2/s -> m2/s
                D[isp][jsp] = Dij;
                D[jsp][isp] = Dij;
            }
        }
    }


private:
    int mNSpecies;
    int mElectronIdx = -1;
    int[] mMolecularSpecies;
    double[] mMolMasses;
    double[] mParticleMass;
    int[] mCharge;
    // working array space
    number[] mMolef;
    number[][] mA_11, mB_11, mC_11, mD_11, mDelta_11, mAlpha;
    number[][] mA_22, mB_22, mC_22, mD_22, mDelta_22, mMu;
    bool[][] mIsCoulombCollision;

    @nogc
    number electronPressureCorrection(ref const(GasState) gs)
    {
    /*
        Compute a correction factor for the collision cross-sections of the charged species,
        based on the electron pressure. This is required as per Gupta, 1989, as explained on
        page 20.

        Notes:
         - Instead of computing p_em, the maximum allowed pressure, we just limit log(Lambda)
           to a minimum of one, which is equivalent, but neater.
         - This function becomes undefined for p_e==0.0, which can happen if no electrons are
           present. Analytically, this would be fine because later on the result will get
           multiplied by the electron mole fraction, which is also equal to zero. But numerically
           the divide be zero causes problems. To fix this we just limit p_e. It should be
           impossible for termTerm or logLambda to overflow as long as p_e is finite, so I chose
           a really small limit.
         - It would be nicer if we could solve this problem by switching to a series expansion
           as p_e gets small, but it messes up the code because you need the expression in terms
           of p_e*log(Lambda).

        @author: Nick Gibbons
    */
        number p_e = fmax(1e-32, gs.p_e);
        number p_e_atm = p_e/P_atm;
        number tempTerm = gs.T/1000.0/pow(p_e_atm, 0.25);
        number logLambda = 0.5*log(2.09e-2*pow(tempTerm, 4.0) + 1.52*pow(tempTerm, 8./3.));
        logLambda = fmax(logLambda, 1.0); // Equivalent to limiting by p_em, from equation 23g
        return logLambda;
    }

    @nogc
    void computeDelta11(ref const(GasState) gs)
    {
        number logLambda = 1.0;
        if (mElectronIdx != -1) logLambda = electronPressureCorrection(gs);

        double kB = Boltzmann_constant;
        number T_CI;
        number log_T_CI;

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                 if (isp != mElectronIdx) {
                    // heavy-particle colliders: use transrotational temperature in calculation
                    T_CI = gs.T;
                    log_T_CI = log(T_CI);
                }
                else {
                    // collisions with electron: use vibroelectronic temperature in calculation
                    T_CI = gs.T_modes[0];
                    log_T_CI = log(T_CI);
                }
                number expnt = mA_11[isp][jsp]*(log_T_CI)^^2 + mB_11[isp][jsp]*log_T_CI + mC_11[isp][jsp];
                number pi_Omega_11 = exp(mD_11[isp][jsp])*pow(T_CI, expnt);
                if (mIsCoulombCollision[isp][jsp])
                    pi_Omega_11 *= logLambda;

                mDelta_11[isp][jsp] = (8.0/3)*1.546e-20*sqrt(2.0*mMu[isp][jsp]/(to!double(PI)*R_universal_cal*T_CI))*pi_Omega_11;
                mDelta_11[jsp][isp] = mDelta_11[isp][jsp];
            }
        }
    }

    @nogc
    void computeDelta22(ref const(GasState) gs)
    {
        number logLambda = 1.0;
        if (mElectronIdx != -1) logLambda = electronPressureCorrection(gs);

        double kB = Boltzmann_constant;
        number T_CI;
        number log_T_CI;

        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. isp+1) {
                if (isp != mElectronIdx) {
                    // heavy-particle colliders: use transrotational temperature in calculation
                    T_CI = gs.T;
                    log_T_CI = log(T_CI);
                }
                else {
                    // collisions with electron: use vibroelectronic temperature in calculation
                    T_CI = gs.T_modes[0];
                    log_T_CI = log(T_CI);
                }
                number expnt = mA_22[isp][jsp]*(log_T_CI)^^2 + mB_22[isp][jsp]*log_T_CI + mC_22[isp][jsp];
                number pi_Omega_22 = exp(mD_22[isp][jsp])*pow(T_CI, expnt);
                if (mIsCoulombCollision[isp][jsp])
                    pi_Omega_22 *= logLambda;

                mDelta_22[isp][jsp] = (16./5)*1.546e-20*sqrt(2.0*mMu[isp][jsp]/(to!double(PI)*R_universal_cal*T_CI))*pi_Omega_22;
                mDelta_22[jsp][isp] = mDelta_22[isp][jsp];
            }
        }
    }
}

version(complex_numbers) {

class CustomLAMPTransProps : TransportPropertiesModel {
public:
    this(lua_State *L, string[] speciesNames)
    {
        throw new Error("CustomLAMP transport_model is not yet available in complex-number LMR builds.");
    }

    void numberDensities(ref const(GasState) gs, ref number[] result)
    {
        throw new Error("customLAMPNumberDensities() is not available in complex-number LMR builds.");
    }

    @nogc
    override void updateTransProps(ref GasState gs, GasModel gm)
    {
        throw new Error("CustomLAMP transport_model is not yet available in complex-number LMR builds.");
    }

    @nogc
    override void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D)
    {
        throw new Error("CustomLAMP transport_model is not yet available in complex-number LMR builds.");
    }
}

} else {

int customLAMPLampSpeciesIndex(string speciesName)
{
    switch (speciesName) {
    case "N2": return 0;
    case "O2": return 1;
    case "N": return 2;
    case "O": return 3;
    case "NO": return 4;
    case "N+": return 5;
    case "O+": return 6;
    case "N2+": return 7;
    case "O2+": return 8;
    case "NO+": return 9;
    case "e-": return 10;
    default: return -1;
    }
}

@nogc double customLAMPLampMolarMass(int lampIndex)
{
    final switch (lampIndex) {
    case 0: return 0.028;       // N2
    case 1: return 0.032;       // O2
    case 2: return 0.014;       // N
    case 3: return 0.016;       // O
    case 4: return 0.030;       // NO
    case 5: return 0.014;       // N+
    case 6: return 0.016;       // O+
    case 7: return 0.028;       // N2+
    case 8: return 0.032;       // O2+
    case 9: return 0.030;       // NO+
    case 10: return 5.484e-07;  // e-
    }
}

@nogc number customLAMPNumberDensity(ref const(GasState) gs, int gdtkIndex, int lampIndex)
{
    enum double numberDensityFloor = 1.0e-10;
    number massf = gs.massf[gdtkIndex];
    if (customLAMPReal(massf) <= 0.0) return customLAMPNumberFromDouble(numberDensityFloor);
    number n = gs.rho * massf / customLAMPLampMolarMass(lampIndex) * CUSTOM_LAMP_AVOGADRO_NUMBER;
    return customLAMPReal(n) > numberDensityFloor ? n : customLAMPNumberFromDouble(numberDensityFloor);
}

@nogc double customLAMPSpeciesMass(int lampIndex)
{
    final switch (lampIndex) {
    case 0: return 28.0 * CUSTOM_LAMP_AMU;
    case 1: return 32.0 * CUSTOM_LAMP_AMU;
    case 2: return 14.0 * CUSTOM_LAMP_AMU;
    case 3: return 16.0 * CUSTOM_LAMP_AMU;
    case 4: return 30.0 * CUSTOM_LAMP_AMU;
    case 5: return 14.0 * CUSTOM_LAMP_AMU;
    case 6: return 16.0 * CUSTOM_LAMP_AMU;
    case 7: return 28.0 * CUSTOM_LAMP_AMU;
    case 8: return 32.0 * CUSTOM_LAMP_AMU;
    case 9: return 30.0 * CUSTOM_LAMP_AMU;
    case 10: return CUSTOM_LAMP_ME;
    }
}

@nogc double customLAMPCharge(int lampIndex)
{
    return lampIndex == 10 ? -1.0 : (lampIndex >= 5 ? 1.0 : 0.0);
}

@nogc double customLAMPLaricchiutaPolarizability(int lampIndex)
{
    final switch (lampIndex) {
    case 0: return 1.749915;
    case 1: return 1.599932;
    case 2: return 1.099759;
    case 3: return 0.799950;
    case 4: return 1.700022;
    case 5: return 0.550850;
    case 6: return 0.279212;
    case 7: return 1.748201;
    case 8: return 0.798385;
    case 9: return 0.851057;
    case 10: return 0.0;
    }
}

@nogc double customLAMPEffectiveElectrons(int lampIndex)
{
    final switch (lampIndex) {
    case 0: return 7.593140;
    case 1: return 9.319243;
    case 2: return 4.355986;
    case 3: return 5.249715;
    case 4: return 8.436012;
    case 5: return 0.0;
    case 6: return 0.0;
    case 7: return 0.0;
    case 8: return 0.0;
    case 9: return 0.0;
    case 10: return 0.0;
    }
}

@nogc double customLAMPScaledIonNeutralPolarizability(int lampIndex)
{
    final switch (lampIndex) {
    case 0: return 1.76;
    case 1: return 1.60;
    case 2: return 1.13;
    case 3: return 0.77;
    case 4: return 1.70;
    case 5: return 0.0;
    case 6: return 0.0;
    case 7: return 0.0;
    case 8: return 0.0;
    case 9: return 0.0;
    case 10: return 0.0;
    }
}

@nogc double customLAMPEqA1(double T, double a1, double a2, double a3, double a4, double a5, double a6)
{
    return (a1 + a2 * pow(T, a3)) / (a4 + a5 * pow(T, a6)) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
}

@nogc double customLAMPEqA2(
    double lnT,
    double a1, double a2, double a3, double a4, double a5,
    double a6, double a7, double a8, double a9)
{
    double exponent = (lnT - a1) / a2;
    double expPos = exp(exponent);
    double expNeg = exp(-exponent);
    return (
        a3 * pow(lnT, a6) * expPos / (expPos + expNeg)
        + a7 * exp(-pow((lnT - a8) / a9, 2.0))
        + a4 * (1.0 + pow(lnT, a5))
    ) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
}

@nogc double customLAMPFermiSwitch(double xLnTStar, double center, double width)
{
    double exponent = (xLnTStar - center) / width;
    double expPos = exp(exponent);
    double expNeg = exp(-exponent);
    return expPos / (expPos + expNeg);
}

@nogc double customLAMPPhenomenologicalA11(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 7.884756e-1 + (-2.438494e-2) * beta;
    case 1: return -2.952759e-1 + (-1.744149e-3) * beta;
    case 2: return 5.020892e-1 + 4.316985e-2 * beta;
    case 3: return -9.042460e-1 + (-4.017103e-2) * beta;
    case 4: return -3.373058e0 + 2.458538e-1 * beta + (-4.850047e-3) * beta * beta;
    case 5: return 4.161981e0 + 2.202737e-1 * beta + (-1.718010e-2) * beta * beta;
    case 6: return 2.462523e0 + 3.231308e-1 * beta + (-2.281072e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalA12(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 7.123565e-1 + (-2.688875e-2) * beta;
    case 1: return -2.910530e-1 + (-2.065175e-3) * beta;
    case 2: return 4.187065e-2 + 4.060236e-2 * beta;
    case 3: return -9.287685e-1 + (-2.342270e-2) * beta;
    case 4: return -3.598542e0 + 2.545120e-1 * beta + (-4.685966e-3) * beta * beta;
    case 5: return 3.934824e0 + 2.699944e-1 * beta + (-2.009886e-2) * beta * beta;
    case 6: return 2.578084e0 + 3.449024e-1 * beta + (-2.292710e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalA13(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 6.606022e-1 + (-2.831448e-2) * beta;
    case 1: return -2.870900e-1 + (-2.232872e-3) * beta;
    case 2: return -2.519690e-1 + 3.778211e-2 * beta;
    case 3: return -9.173046e-1 + (-1.864476e-2) * beta;
    case 4: return -3.776812e0 + 2.552528e-1 * beta + (-4.237220e-3) * beta * beta;
    case 5: return 3.768103e0 + 3.155025e-1 * beta + (-2.218849e-2) * beta * beta;
    case 6: return 2.695440e0 + 3.597998e-1 * beta + (-2.267102e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalA22(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 7.898524e-1 + (-2.114115e-2) * beta;
    case 1: return -2.998325e-1 + (-1.243977e-3) * beta;
    case 2: return 7.077103e-1 + 3.583907e-2 * beta;
    case 3: return -8.946857e-1 + (-2.473947e-2) * beta;
    case 4: return -2.958969e0 + 2.303358e-1 * beta + (-5.226562e-3) * beta * beta;
    case 5: return 4.348412e0 + 1.920321e-1 * beta + (-1.496557e-2) * beta * beta;
    case 6: return 2.205440e0 + 2.567027e-1 * beta + (-1.861359e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalNI_A11(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 9.851755e-1 + (-2.870704e-2) * beta;
    case 1: return -4.737800e-1 + (-1.370344e-3) * beta;
    case 2: return 7.080799e-1 + 4.575312e-3 * beta;
    case 3: return -1.239439e0 + (-3.683605e-2) * beta;
    case 4: return -4.638467e0 + 4.418904e-1 * beta + (-1.220292e-2) * beta * beta;
    case 5: return 3.841835e0 + 3.277658e-1 * beta + (-2.660275e-2) * beta * beta;
    case 6: return 2.317342e0 + 3.912768e-1 * beta + (-3.136223e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalNI_A12(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 8.361751e-1 + (-3.201292e-2) * beta;
    case 1: return -4.707355e-1 + (-1.783284e-3) * beta;
    case 2: return 1.771157e-1 + 1.172773e-2 * beta;
    case 3: return -1.094937e0 + (-3.115598e-2) * beta;
    case 4: return -4.976384e0 + 4.708074e-1 * beta + (-1.283818e-2) * beta * beta;
    case 5: return 3.645873e0 + 3.699452e-1 * beta + (-2.988684e-2) * beta * beta;
    case 6: return 2.428864e0 + 4.267351e-1 * beta + (-3.278874e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalNI_A13(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 7.440562e-1 + (-3.453851e-2) * beta;
    case 1: return -4.656306e-1 + (-2.097901e-3) * beta;
    case 2: return -1.465752e-1 + 1.446209e-2 * beta;
    case 3: return -1.080410e0 + (-2.712029e-2) * beta;
    case 4: return -5.233907e0 + 4.846691e-1 * beta + (-1.280346e-2) * beta * beta;
    case 5: return 3.489814e0 + 4.140270e-1 * beta + (-3.250138e-2) * beta * beta;
    case 6: return 2.529678e0 + 4.515088e-1 * beta + (-3.339293e-2) * beta * beta;
    }
}

@nogc double customLAMPPhenomenologicalNI_A22(int rowIndex, double beta)
{
    final switch (rowIndex) {
    case 0: return 9.124518e-1 + (-2.398461e-2) * beta;
    case 1: return -4.697184e-1 + (-7.809681e-4) * beta;
    case 2: return 1.031053e0 + 4.069668e-3 * beta;
    case 3: return -1.090782e0 + (-2.413508e-2) * beta;
    case 4: return -4.127243e0 + 4.302667e-1 * beta + (-1.352874e-2) * beta * beta;
    case 5: return 4.059078e0 + 2.597379e-1 * beta + (-2.169951e-2) * beta * beta;
    case 6: return 2.086906e0 + 2.920310e-1 * beta + (-2.560437e-2) * beta * beta;
    }
}

@nogc double customLAMPLnOmegaStar(double x, double a1, double a2, double a3, double a4, double a5, double a6, double a7)
{
    return (a1 + a2 * x) * customLAMPFermiSwitch(x, a3, a4)
        + a5 * customLAMPFermiSwitch(x, a6, a7);
}

@nogc double customLAMPPhenomenologicalNN(double T, int speciesI, int speciesJ)
{
    double alphaI = customLAMPLaricchiutaPolarizability(speciesI);
    double alphaJ = customLAMPLaricchiutaPolarizability(speciesJ);
    double neffI = customLAMPEffectiveElectrons(speciesI);
    double neffJ = customLAMPEffectiveElectrons(speciesJ);
    double beta = 6.0 + 5.0 / (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0));
    double re = 1.767 * (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0))
        / pow(alphaI * alphaJ, 9.5e-2);
    double cd = 1.57e1 * alphaI * alphaJ / (sqrt(alphaI / neffI) + sqrt(alphaJ / neffJ));
    double eps0 = 7.2e-1 * cd / pow(re, 6.0);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 8.002000e-1 * pow(beta, 4.925600e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    double a1 = customLAMPPhenomenologicalA11(0, beta);
    double a2 = customLAMPPhenomenologicalA11(1, beta);
    double a3 = customLAMPPhenomenologicalA11(2, beta);
    double a4 = customLAMPPhenomenologicalA11(3, beta);
    double a5 = customLAMPPhenomenologicalA11(4, beta);
    double a6 = customLAMPPhenomenologicalA11(5, beta);
    double a7 = customLAMPPhenomenologicalA11(6, beta);
    double lnOmegaStar = customLAMPLnOmegaStar(x, a1, a2, a3, a4, a5, a6, a7);
    return exp(lnOmegaStar) * scale;
}

@nogc double customLAMPPhenomenologicalNN22(double T, int speciesI, int speciesJ)
{
    double alphaI = customLAMPLaricchiutaPolarizability(speciesI);
    double alphaJ = customLAMPLaricchiutaPolarizability(speciesJ);
    double neffI = customLAMPEffectiveElectrons(speciesI);
    double neffJ = customLAMPEffectiveElectrons(speciesJ);
    double beta = 6.0 + 5.0 / (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0));
    double re = 1.767 * (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0))
        / pow(alphaI * alphaJ, 9.5e-2);
    double cd = 1.57e1 * alphaI * alphaJ / (sqrt(alphaI / neffI) + sqrt(alphaJ / neffJ));
    double eps0 = 7.2e-1 * cd / pow(re, 6.0);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 8.002000e-1 * pow(beta, 4.925600e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    double a1 = customLAMPPhenomenologicalA22(0, beta);
    double a2 = customLAMPPhenomenologicalA22(1, beta);
    double a3 = customLAMPPhenomenologicalA22(2, beta);
    double a4 = customLAMPPhenomenologicalA22(3, beta);
    double a5 = customLAMPPhenomenologicalA22(4, beta);
    double a6 = customLAMPPhenomenologicalA22(5, beta);
    double a7 = customLAMPPhenomenologicalA22(6, beta);
    return exp(customLAMPLnOmegaStar(x, a1, a2, a3, a4, a5, a6, a7)) * scale;
}

@nogc double customLAMPPhenomenologicalNN12(double T, int speciesI, int speciesJ)
{
    double alphaI = customLAMPLaricchiutaPolarizability(speciesI);
    double alphaJ = customLAMPLaricchiutaPolarizability(speciesJ);
    double neffI = customLAMPEffectiveElectrons(speciesI);
    double neffJ = customLAMPEffectiveElectrons(speciesJ);
    double beta = 6.0 + 5.0 / (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0));
    double re = 1.767 * (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0))
        / pow(alphaI * alphaJ, 9.5e-2);
    double cd = 1.57e1 * alphaI * alphaJ / (sqrt(alphaI / neffI) + sqrt(alphaJ / neffJ));
    double eps0 = 7.2e-1 * cd / pow(re, 6.0);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 8.002000e-1 * pow(beta, 4.925600e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    return exp(customLAMPLnOmegaStar(
        x,
        customLAMPPhenomenologicalA12(0, beta),
        customLAMPPhenomenologicalA12(1, beta),
        customLAMPPhenomenologicalA12(2, beta),
        customLAMPPhenomenologicalA12(3, beta),
        customLAMPPhenomenologicalA12(4, beta),
        customLAMPPhenomenologicalA12(5, beta),
        customLAMPPhenomenologicalA12(6, beta))) * scale;
}

@nogc double customLAMPPhenomenologicalNN13(double T, int speciesI, int speciesJ)
{
    double alphaI = customLAMPLaricchiutaPolarizability(speciesI);
    double alphaJ = customLAMPLaricchiutaPolarizability(speciesJ);
    double neffI = customLAMPEffectiveElectrons(speciesI);
    double neffJ = customLAMPEffectiveElectrons(speciesJ);
    double beta = 6.0 + 5.0 / (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0));
    double re = 1.767 * (pow(alphaI, 1.0 / 3.0) + pow(alphaJ, 1.0 / 3.0))
        / pow(alphaI * alphaJ, 9.5e-2);
    double cd = 1.57e1 * alphaI * alphaJ / (sqrt(alphaI / neffI) + sqrt(alphaJ / neffJ));
    double eps0 = 7.2e-1 * cd / pow(re, 6.0);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 8.002000e-1 * pow(beta, 4.925600e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    return exp(customLAMPLnOmegaStar(
        x,
        customLAMPPhenomenologicalA13(0, beta),
        customLAMPPhenomenologicalA13(1, beta),
        customLAMPPhenomenologicalA13(2, beta),
        customLAMPPhenomenologicalA13(3, beta),
        customLAMPPhenomenologicalA13(4, beta),
        customLAMPPhenomenologicalA13(5, beta),
        customLAMPPhenomenologicalA13(6, beta))) * scale;
}

@nogc double customLAMPPhenomenologicalNI(double T, int neutralSpecies, int ionSpecies, bool omega22)
{
    double alphaNeutral = customLAMPLaricchiutaPolarizability(neutralSpecies);
    double alphaIon = customLAMPLaricchiutaPolarizability(ionSpecies);
    double beta = 6.0 + 5.0 / (pow(alphaNeutral, 1.0 / 3.0) + pow(alphaIon, 1.0 / 3.0));
    double rho = alphaIon / ((1.0 + pow(2.0 * alphaIon / alphaNeutral, 2.0 / 3.0)) * sqrt(alphaNeutral));
    double re = 1.767 * (pow(alphaNeutral, 1.0 / 3.0) + pow(alphaIon, 1.0 / 3.0))
        / pow(alphaNeutral * alphaIon * (1.0 + 1.0 / rho), 9.5e-2);
    double eps0 = 5.2e0 * alphaNeutral / pow(re, 4.0) * (1.0 + rho);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 7.564000e-1 * pow(beta, 6.460500e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    double a1 = omega22 ? customLAMPPhenomenologicalNI_A22(0, beta) : customLAMPPhenomenologicalNI_A11(0, beta);
    double a2 = omega22 ? customLAMPPhenomenologicalNI_A22(1, beta) : customLAMPPhenomenologicalNI_A11(1, beta);
    double a3 = omega22 ? customLAMPPhenomenologicalNI_A22(2, beta) : customLAMPPhenomenologicalNI_A11(2, beta);
    double a4 = omega22 ? customLAMPPhenomenologicalNI_A22(3, beta) : customLAMPPhenomenologicalNI_A11(3, beta);
    double a5 = omega22 ? customLAMPPhenomenologicalNI_A22(4, beta) : customLAMPPhenomenologicalNI_A11(4, beta);
    double a6 = omega22 ? customLAMPPhenomenologicalNI_A22(5, beta) : customLAMPPhenomenologicalNI_A11(5, beta);
    double a7 = omega22 ? customLAMPPhenomenologicalNI_A22(6, beta) : customLAMPPhenomenologicalNI_A11(6, beta);
    return exp(customLAMPLnOmegaStar(x, a1, a2, a3, a4, a5, a6, a7)) * scale;
}

@nogc double customLAMPPhenomenologicalNIByOrder(double T, int neutralSpecies, int ionSpecies, int order)
{
    double alphaNeutral = customLAMPLaricchiutaPolarizability(neutralSpecies);
    double alphaIon = customLAMPLaricchiutaPolarizability(ionSpecies);
    double beta = 6.0 + 5.0 / (pow(alphaNeutral, 1.0 / 3.0) + pow(alphaIon, 1.0 / 3.0));
    double rho = alphaIon / ((1.0 + pow(2.0 * alphaIon / alphaNeutral, 2.0 / 3.0)) * sqrt(alphaNeutral));
    double re = 1.767 * (pow(alphaNeutral, 1.0 / 3.0) + pow(alphaIon, 1.0 / 3.0))
        / pow(alphaNeutral * alphaIon * (1.0 + 1.0 / rho), 9.5e-2);
    double eps0 = 5.2e0 * alphaNeutral / pow(re, 4.0) * (1.0 + rho);
    double x = log(CUSTOM_LAMP_KB * T / (eps0 * CUSTOM_LAMP_EC));
    double x0 = 7.564000e-1 * pow(beta, 6.460500e-2);
    double scale = pow(x0 * re, 2.0) * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
    double a1, a2, a3, a4, a5, a6, a7;
    if (order == 12) {
        a1 = customLAMPPhenomenologicalNI_A12(0, beta);
        a2 = customLAMPPhenomenologicalNI_A12(1, beta);
        a3 = customLAMPPhenomenologicalNI_A12(2, beta);
        a4 = customLAMPPhenomenologicalNI_A12(3, beta);
        a5 = customLAMPPhenomenologicalNI_A12(4, beta);
        a6 = customLAMPPhenomenologicalNI_A12(5, beta);
        a7 = customLAMPPhenomenologicalNI_A12(6, beta);
    } else {
        a1 = customLAMPPhenomenologicalNI_A13(0, beta);
        a2 = customLAMPPhenomenologicalNI_A13(1, beta);
        a3 = customLAMPPhenomenologicalNI_A13(2, beta);
        a4 = customLAMPPhenomenologicalNI_A13(3, beta);
        a5 = customLAMPPhenomenologicalNI_A13(4, beta);
        a6 = customLAMPPhenomenologicalNI_A13(5, beta);
        a7 = customLAMPPhenomenologicalNI_A13(6, beta);
    }
    return exp(customLAMPLnOmegaStar(x, a1, a2, a3, a4, a5, a6, a7)) * scale;
}

@nogc double customLAMPScaledPairOmega11(double T, int neutralSpeciesIndex)
{
    return 425.4 * sqrt(customLAMPScaledIonNeutralPolarizability(neutralSpeciesIndex) / T)
        * CUSTOM_LAMP_CS * CUSTOM_LAMP_PI;
}

@nogc double customLAMPScaledPairOmega22(double T, int neutralSpeciesIndex)
{
    return 0.8710 * customLAMPScaledPairOmega11(T, neutralSpeciesIndex);
}

@nogc double customLAMPScaledPairOmega12(double T, int neutralSpeciesIndex)
{
    return 0.8333 * customLAMPScaledPairOmega11(T, neutralSpeciesIndex);
}

@nogc double customLAMPScaledPairOmega13(double T, int neutralSpeciesIndex)
{
    return 0.7292 * customLAMPScaledPairOmega11(T, neutralSpeciesIndex);
}

@nogc double customLAMPElectronNeutralOmega(
    double Te,
    int neutralSpeciesIndex,
    int order)
{
    double lnTe = log(Te);
    if (neutralSpeciesIndex == 0) {
        if (order == 11) return customLAMPEqA2(lnTe, 4.4763e0, 8.9123e0, 8.6558e1, -1.1907e1, 7.3977e-1, 8.1569e-2, 1.7704e0, 9.2013e0, 8.0978e-1);
        if (order == 12) return customLAMPEqA2(lnTe, 6.2997e0, 2.8143e0, 6.0373e0, -6.1715e-7, 5.9733e0, -1.3242e-1, 2.1061e0, 8.9217e0, 7.2869e-1);
        if (order == 13) return customLAMPEqA2(lnTe, 6.0853e0, 2.5330e0, 9.8515e0, -1.1814e-10, 9.2906e0, -4.0522e-1, 2.4046e0, 8.6973e0, 6.7041e-1);
        if (order == 14) return customLAMPEqA2(lnTe, 5.8136e0, 2.4462e0, 9.3499e0, -1.1687e-11, 1.0256e1, -3.9827e-1, 2.6488e0, 8.5107e0, 6.2512e-1);
        if (order == 15) return customLAMPEqA2(lnTe, 5.4329e0, 2.3732e0, 6.5229e0, -5.6704e-11, 9.6931e0, -2.4519e-1, 2.8589e0, 8.3525e0, 5.9211e-1);
        if (order == 22) return customLAMPEqA2(lnTe, 7.7099e0, 3.6418e0, 6.9892e1, -2.9740e-1, 1.3600e0, -7.0673e-1, 2.1153e0, 8.9418e0, 8.3996e-1);
    } else if (neutralSpeciesIndex == 1) {
        if (order == 11) return customLAMPEqA2(lnTe, 8.3307e0, 2.1057e0, 1.3703e2, 1.6462e-3, 2.2034e0, -1.6638e0, -7.1918e-1, 9.7316e0, 1.0537e0);
        if (order == 12) return customLAMPEqA2(lnTe, 8.3440e0, 2.1878e0, 2.6005e2, -4.7740e-2, -7.8236e-1, -1.8609e0, -9.0948e-1, 9.4672e0, 1.0447e0);
        if (order == 13) return customLAMPEqA2(lnTe, 8.2363e0, 2.2106e0, 3.5049e2, -1.0523e-1, -8.9547e-1, -1.9861e0, -9.7555e-1, 9.2465e0, 1.0010e0);
        if (order == 14) return customLAMPEqA2(lnTe, 7.7713e0, 2.2180e0, 8.0242e1, -2.9025e-6, 4.7426e0, -1.3707e0, -9.8561e-1, 9.0635e0, 9.5788e-1);
        if (order == 15) return customLAMPEqA2(lnTe, 8.0580e0, 2.2437e0, 4.9008e2, -2.0701e-1, -1.0225e0, -2.1302e0, -1.0577e0, 8.9127e0, 9.4543e-1);
        if (order == 22) return customLAMPEqA2(lnTe, 8.1575e0, 2.8053e0, 9.3652e0, -4.1748e-4, 3.6655e0, -1.8979e-1, -7.6150e-1, 9.3932e0, 9.6517e-1);
    } else if (neutralSpeciesIndex == 2) {
        if (order == 11) return customLAMPEqA2(lnTe, 8.9279e0, 1.1613e0, 1.6494e0, 1.2796e-1, 2.7394e-1, -2.7594e-1, 0.0000e0, 0.0000e0, 1.0000e0);
        if (order == 12) return customLAMPEqA2(lnTe, 8.6314e0, 1.0738e0, 1.8882e0, 1.2214e-1, 3.3398e-1, -3.4447e-1, 0.0000e0, 0.0000e0, 1.0000e0);
        if (order == 13) return customLAMPEqA2(lnTe, 8.4009e0, 9.8544e-1, 1.9023e0, 1.1026e-1, 4.4086e-1, -3.6762e-1, 0.0000e0, 0.0000e0, 1.0000e0);
        if (order == 14) return customLAMPEqA2(lnTe, 8.2043e0, 9.4398e-1, 1.8527e0, 1.0647e-1, 4.8117e-1, -3.6661e-1, 0.0000e0, 0.0000e0, 1.0000e0);
        if (order == 15) return customLAMPEqA2(lnTe, 8.0388e0, 9.0527e-1, 1.8833e0, 1.0051e-1, 5.4116e-1, -3.8773e-1, 0.0000e0, 0.0000e0, 1.0000e0);
        if (order == 22) return customLAMPEqA2(lnTe, 8.2257e0, 9.5133e-1, 8.7677e-1, 6.8519e-2, 1.6570e-1, 3.4127e-2, 3.5900e-1, 9.6380e0, 9.0876e-1);
    } else if (neutralSpeciesIndex == 3) {
        if (order == 11) return customLAMPEqA2(lnTe, 1.0568e1, 1.3640e0, 3.0548e4, 2.4439e-1, 1.2512e-1, -4.2290e0, 2.7943e-1, 8.6407e0, 1.7624e0);
        if (order == 12) return customLAMPEqA2(lnTe, 1.0276e1, 1.3728e0, 7.8394e4, 2.2579e-1, 2.0701e-1, -4.6626e0, 2.4241e-1, 8.1977e0, 1.5622e0);
        if (order == 13) return customLAMPEqA2(lnTe, 1.0069e1, 1.2781e0, 1.2008e5, 2.1584e-1, 2.7462e-1, -4.9024e0, 2.5935e-1, 8.0276e0, 1.5527e0);
        if (order == 14) return customLAMPEqA2(lnTe, 9.9739e0, 1.3564e0, 2.7669e5, 2.2906e-1, 1.7466e-1, -5.2350e0, 2.4244e-1, 7.7435e0, 1.4870e0);
        if (order == 15) return customLAMPEqA2(lnTe, 9.9629e0, 1.4615e0, 9.4493e5, 1.9454e-1, 2.2494e-1, -5.7280e0, 1.9924e-1, 7.4822e0, 1.3722e0);
        if (order == 22) return customLAMPEqA2(lnTe, 9.9767e0, 1.6871e0, 2.9697e5, 5.6206e-2, 1.4522e-1, -5.1321e0, 7.1705e-2, 8.4309e0, 1.3830e0);
    } else if (neutralSpeciesIndex == 4) {
        if (order == 11) return customLAMPEqA2(lnTe, -7.2965e-1, 1.3332e1, 7.8565e2, -1.9259e2, 2.8336e-1, -6.3703e-2, 9.3029e-1, 8.8086e0, 9.6216e-1);
        if (order == 12) return customLAMPEqA2(lnTe, -1.4887e1, 1.8148e1, 2.5286e3, -1.0586e3, 2.5251e-1, 9.4792e-2, 1.0562e0, 8.5064e0, 8.7879e-1);
        if (order == 13) return customLAMPEqA2(lnTe, -9.5704e0, 1.7653e1, 1.4649e3, -5.5148e2, 3.6509e-1, 1.4218e-1, 1.1651e0, 8.2614e0, 8.2346e-1);
        if (order == 14) return customLAMPEqA2(lnTe, -4.8381e0, 1.6713e1, 8.8349e2, -2.8954e2, 5.0399e-1, 2.0687e-1, 1.2667e0, 8.0651e0, 7.9767e-1);
        if (order == 15) return customLAMPEqA2(lnTe, 5.8041e-1, 1.3987e1, 3.9144e2, -9.7881e1, 7.0323e-1, 2.8223e-1, 1.3560e0, 7.9087e0, 7.8024e-1);
        if (order == 22) return customLAMPEqA2(lnTe, 3.8989e-2, 1.4692e1, 6.6332e2, -1.6621e2, 5.9703e-1, 1.9528e-1, 1.0178e0, 8.6724e0, 9.6611e-1);
    }
    return 0.0;
}

@nogc double customLAMPElectronNeutralOmega11(double Te, int neutralSpeciesIndex)
{
    return customLAMPElectronNeutralOmega(Te, neutralSpeciesIndex, 11);
}

@nogc double customLAMPEmpiricalNeutralOmega11(double T, int speciesI, int speciesJ)
{
    int i = speciesI < speciesJ ? speciesI : speciesJ;
    int j = speciesI < speciesJ ? speciesJ : speciesI;
    if (i == 0 && j == 4) return customLAMPEqA1(T, 1.6324e2, 7.2640e3, -9.5137e-1, 8.2746e0, 1.0000e0, 3.4180e-1);
    if (i == 1 && j == 4) return customLAMPEqA1(T, 1.5090e2, 6.7215e3, -9.1171e-1, 9.6404e0, 8.2038e-1, 3.4903e-1);
    if (i == 2 && j == 2) return customLAMPEqA1(T, 1.9939e1, 9.0652e2, -8.4875e-1, 2.3809e0, 3.8782e-2, 4.8674e-1);
    if (i == 2 && j == 3) return customLAMPEqA1(T, 2.2683e1, 1.2470e3, -1.0271e0, 1.9320e0, 6.4407e-2, 5.0887e-1);
    if (i == 2 && j == 4) return customLAMPEqA1(T, 1.8376e1, 5.8513e2, -1.0885e0, 2.1989e-2, 3.7017e-1, 2.7990e-1);
    if (i == 3 && j == 4) return customLAMPEqA1(T, 1.0948e2, 3.9522e3, -9.2591e-1, 5.6774e0, 1.0000e0, 3.4478e-1);
    if (i == 4 && j == 4) return customLAMPEqA1(T, 3.5215e2, 1.7562e4, -8.8485e-1, 2.9253e1, 6.3570e-1, 4.5570e-1);
    throw new Error("unsupported CustomLAMP empirical neutral pair");
}

@nogc double customLAMPEmpiricalNeutralOmega22(double T, int speciesI, int speciesJ)
{
    int i = speciesI < speciesJ ? speciesI : speciesJ;
    int j = speciesI < speciesJ ? speciesJ : speciesI;
    if (i == 0 && j == 4) return customLAMPEqA1(T, 1.7855e2, 7.4922e3, -8.3567e-1, 1.3283e1, 4.0509e-1, 3.9878e-1);
    if (i == 1 && j == 4) return customLAMPEqA1(T, 2.4059e2, 1.0381e4, -8.1156e-1, 2.1388e1, 3.9864e-1, 4.2191e-1);
    if (i == 2 && j == 2) return customLAMPEqA1(T, 2.3736e1, 9.1364e2, -7.7981e-1, 2.7874e0, 4.9002e-2, 4.6224e-1);
    if (i == 2 && j == 3) return customLAMPEqA1(T, 2.1163e1, 8.6509e2, -8.1169e-1, 2.6270e0, 2.2254e-2, 5.6508e-1);
    if (i == 2 && j == 4) return customLAMPEqA1(T, 5.2091e0, 1.2579e2, -7.9721e-1, 2.2715e-1, 6.0530e-2, 2.9981e-1);
    if (i == 3 && j == 4) return customLAMPEqA1(T, 3.4111e2, 1.2306e4, -7.6446e-1, 3.4503e1, 9.6152e-1, 4.1455e-1);
    if (i == 4 && j == 4) return customLAMPEqA1(T, 3.7951e2, 1.7383e4, -7.9157e-1, 3.7919e1, 1.8054e-1, 5.4195e-1);
    throw new Error("unsupported CustomLAMP empirical neutral pair");
}

@nogc double customLAMPEmpiricalNeutralOmega12(double T, int speciesI, int speciesJ)
{
    int i = speciesI < speciesJ ? speciesI : speciesJ;
    int j = speciesI < speciesJ ? speciesJ : speciesI;
    if (i == 0 && j == 4) return customLAMPEqA1(T, 1.1646e2, 1.9357e4, -1.4606e0, 4.4626e0, 8.8277e-1, 3.3773e-1);
    if (i == 1 && j == 4) return customLAMPEqA1(T, 9.6359e1, 1.1890e4, -1.3477e0, 4.4156e0, 7.4636e-1, 3.3393e-1);
    if (i == 2 && j == 2) return customLAMPEqA1(T, 9.3402e0, 1.0845e3, -1.2761e0, 7.9736e-1, 3.5324e-2, 4.4711e-1);
    if (i == 2 && j == 3) return customLAMPEqA1(T, 1.0660e1, 2.3223e3, -1.4924e0, 9.1182e-1, 2.4379e-2, 5.4850e-1);
    if (i == 2 && j == 4) return customLAMPEqA1(T, 1.9001e0, 6.2564e2, -1.7235e0, 4.7035e-2, 2.4421e-2, 3.2671e-1);
    if (i == 3 && j == 4) return customLAMPEqA1(T, 4.9973e1, 6.8422e3, -1.4022e0, 2.4320e0, 4.1518e-1, 3.6638e-1);
    if (i == 4 && j == 4) return customLAMPEqA1(T, 1.1494e2, 1.2389e4, -1.2491e0, 7.3328e0, 3.4827e-1, 4.2860e-1);
    throw new Error("unsupported CustomLAMP empirical neutral pair");
}

@nogc double customLAMPEmpiricalNeutralOmega13(double T, int speciesI, int speciesJ)
{
    int i = speciesI < speciesJ ? speciesI : speciesJ;
    int j = speciesI < speciesJ ? speciesJ : speciesI;
    if (i == 0 && j == 4) return customLAMPEqA1(T, 2.6945e1, 1.7648e4, -1.9074e0, 1.1168e0, 1.8471e-1, 3.5476e-1);
    if (i == 1 && j == 4) return customLAMPEqA1(T, 3.2746e1, 1.4087e4, -1.7751e0, 1.4877e0, 2.5246e-1, 3.4294e-1);
    if (i == 2 && j == 2) return customLAMPEqA1(T, 5.5929e0, 2.1522e3, -1.6966e0, 4.4782e-1, 2.4128e-2, 4.4567e-1);
    if (i == 2 && j == 3) return customLAMPEqA1(T, 5.8682e1, 5.0777e4, -1.9247e0, 5.3712e0, 1.0570e-1, 5.8600e-1);
    if (i == 2 && j == 4) return customLAMPEqA1(T, 6.3924e-1, 1.7922e3, -2.3300e0, 2.4350e-2, 6.2276e-3, 3.5811e-1);
    if (i == 3 && j == 4) return customLAMPEqA1(T, 1.5346e1, 1.2856e4, -1.9585e0, 8.1792e-1, 1.1084e-1, 3.8800e-1);
    if (i == 4 && j == 4) return customLAMPEqA1(T, 1.9167e2, 5.3633e4, -1.6030e0, 1.1561e1, 6.3881e-1, 4.3160e-1);
    throw new Error("unsupported CustomLAMP empirical neutral pair");
}

@nogc double customLAMPEmpiricalIonNeutralOmega11(double T, int neutralSpecies, int ionSpecies)
{
    if (neutralSpecies == 0 && ionSpecies == 7) return customLAMPEqA1(T, -1.5707e2, 1.6358e2, -1.5782e-3, 4.6818e-2, 3.7122e-2, 5.8367e-2);
    if (neutralSpecies == 1 && ionSpecies == 8) return customLAMPEqA1(T, -1.1575e2, 3.7386e2, -4.8436e-2, 3.6486e0, 2.6117e-1, 1.4145e-1);
    if (neutralSpecies == 2 && ionSpecies == 5) return customLAMPEqA1(T, -1.3960e2, 3.6208e2, -4.3674e-2, 2.7103e0, 3.0103e-1, 1.2418e-1);
    if (neutralSpecies == 3 && ionSpecies == 6) return customLAMPEqA1(T, -1.1682e2, 3.2563e2, -3.9601e-2, 4.1463e0, 3.1988e-1, 1.0766e-1);
    if (neutralSpecies == 4 && ionSpecies == 9) return customLAMPEqA1(T, -8.1007e2, 8.1015e2, -4.6391e-6, 7.5318e-4, 5.5734e-4, 7.3176e-2);
    throw new Error("unsupported CustomLAMP empirical neutral-ion pair");
}

@nogc double customLAMPEmpiricalIonNeutralOmega12(double T, int neutralSpecies, int ionSpecies)
{
    if (neutralSpecies == 0 && ionSpecies == 7) return customLAMPEqA1(T, -1.5988e2, 1.6168e2, -4.4121e-4, 1.2796e-2, 1.0896e-2, 5.7820e-2);
    if (neutralSpecies == 1 && ionSpecies == 8) return customLAMPEqA1(T, -2.5560e2, 2.7178e2, -2.8455e-3, 1.4886e-1, 1.0396e-1, 7.4847e-2);
    if (neutralSpecies == 2 && ionSpecies == 5) return customLAMPEqA1(T, -1.2812e2, 3.8252e2, -5.0422e-2, 3.2848e0, 2.4802e-1, 1.4072e-1);
    if (neutralSpecies == 3 && ionSpecies == 6) return customLAMPEqA1(T, -1.2217e2, 3.2403e2, -3.6500e-2, 3.9708e0, 4.5589e-1, 1.0641e-1);
    if (neutralSpecies == 4 && ionSpecies == 9) return customLAMPEqA1(T, -8.1015e2, 8.1023e2, -4.9170e-6, 7.6379e-4, 6.2676e-4, 6.9886e-2);
    throw new Error("unsupported CustomLAMP empirical neutral-ion pair");
}

@nogc double customLAMPEmpiricalIonNeutralOmega13(double T, int neutralSpecies, int ionSpecies)
{
    if (neutralSpecies == 0 && ionSpecies == 7) return customLAMPEqA1(T, -1.5655e2, 1.6306e2, -1.6250e-3, 4.7258e-2, 3.9996e-2, 5.7329e-2);
    if (neutralSpecies == 1 && ionSpecies == 8) return customLAMPEqA1(T, -2.0496e2, 5.3410e2, -4.3204e-2, 4.7254e0, 5.0113e-1, 1.1552e-1);
    if (neutralSpecies == 2 && ionSpecies == 5) return customLAMPEqA1(T, -1.3692e2, 3.7076e2, -4.6642e-2, 3.0216e0, 3.0143e-1, 1.3162e-1);
    if (neutralSpecies == 3 && ionSpecies == 6) return customLAMPEqA1(T, -1.1681e1, 2.7424e2, -6.0360e-2, 5.5075e0, 3.2206e-1, 1.7246e-1);
    if (neutralSpecies == 4 && ionSpecies == 9) return customLAMPEqA1(T, -8.1009e2, 8.1017e2, -4.7522e-6, 7.7546e-4, 5.9147e-4, 7.3438e-2);
    throw new Error("unsupported CustomLAMP empirical neutral-ion pair");
}

@nogc double customLAMPEmpiricalIonNeutralOmega22(double T, int neutralSpecies, int ionSpecies)
{
    if (neutralSpecies == 0 && ionSpecies == 7) return customLAMPEqA1(T, -6.1066e1, 1.0830e2, 1.6620e-5, 1.9752e-4, 9.6025e-2, 5.0012e-1);
    if (neutralSpecies == 1 && ionSpecies == 8) return customLAMPEqA1(T, 1.4651e3, 4.5141e2, -4.7817e-1, 1.0000e0, 3.1366e0, 4.9971e-1);
    if (neutralSpecies == 2 && ionSpecies == 5) return customLAMPEqA1(T, -1.85975093e8, 1.85975492e8, 6.4006e-8, 1.6400e1, 1.0000e0, 4.7003e-1);
    if (neutralSpecies == 3 && ionSpecies == 6) return customLAMPEqA1(T, 3.6880e-12, 9.1735e-12, -2.5404e-1, 4.1706e-13, 4.2236e-15, 5.4721e-1);
    if (neutralSpecies == 4 && ionSpecies == 9) return customLAMPEqA1(T, -8.0815e2, 8.0826e2, 1.4346e-9, 3.3750e-7, 2.1668e-4, 5.0007e-1);
    throw new Error("unsupported CustomLAMP empirical neutral-ion pair");
}

@nogc double customLAMPNONOpOmega11(double T)
{
    return customLAMPEqA1(T, -8.1007e2, 8.1015e2, -4.6391e-6, 7.5318e-4, 5.5734e-4, 7.3176e-2);
}

@nogc double customLAMPCoulombScreenedLogTerm(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ,
    out double rhoPerp)
{
    double massI = customLAMPSpeciesMass(lampI);
    double massJ = customLAMPSpeciesMass(lampJ);
    double reducedMass = massI * massJ / (massI + massJ);
    double debyeLength = sqrt(CUSTOM_LAMP_EPS0 * CUSTOM_LAMP_KB * debyeTemperature
        / (CUSTOM_LAMP_EC * CUSTOM_LAMP_EC * electronDensity));
    double thermalRelativeSpeedSq = 2.0 * CUSTOM_LAMP_KB * screeningTemperature / reducedMass;
    double thermalRelativeSpeed = sqrt(thermalRelativeSpeedSq);
    double lambdaDeBroglie = CUSTOM_LAMP_H / (4.0 * CUSTOM_LAMP_PI * reducedMass * thermalRelativeSpeed);
    rhoPerp = abs(customLAMPCharge(lampI) * customLAMPCharge(lampJ))
        * CUSTOM_LAMP_EC * CUSTOM_LAMP_EC
        / (4.0 * CUSTOM_LAMP_PI * CUSTOM_LAMP_EPS0 * reducedMass * thermalRelativeSpeedSq);
    double ionSphereRadius = pow(3.0 / (4.0 * CUSTOM_LAMP_PI * electronDensity), 1.0 / 3.0);
    double bMin = sqrt(lambdaDeBroglie * lambdaDeBroglie + rhoPerp * rhoPerp);
    double bMax = sqrt(debyeLength * debyeLength + ionSphereRadius * ionSphereRadius);
    return 0.5 * log(1.0 + (bMax * bMax) / (bMin * bMin))
        + log(2.0) - 2.0 * CUSTOM_LAMP_EULER_GAMMA;
}

@nogc double customLAMPCoulombOmega11(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega1sFromLog(1, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombPsi(int order)
{
    double result = 0.0;
    foreach (index; 1 .. order) result += 1.0 / to!double(index);
    return result;
}

@nogc double customLAMPCoulombOmega1sFromLog(int order, double screenedLogTerm, double rhoPerp)
{
    return CUSTOM_LAMP_PI * 4.0 / (to!double(order) * to!double(order + 1))
        * rhoPerp * rhoPerp * (screenedLogTerm + customLAMPCoulombPsi(order) - 0.5);
}

@nogc double customLAMPCoulombOmega12(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega1sFromLog(2, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombOmega13(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega1sFromLog(3, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombOmega14(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega1sFromLog(4, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombOmega15(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega1sFromLog(5, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombOmega2sFromLog(int order, double screenedLogTerm, double rhoPerp)
{
    return CUSTOM_LAMP_PI * 12.0 / (to!double(order) * to!double(order + 1))
        * rhoPerp * rhoPerp * (screenedLogTerm + customLAMPCoulombPsi(order) - 1.0);
}

@nogc double customLAMPCoulombOmega22(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return CUSTOM_LAMP_PI * 2.0 * rhoPerp * rhoPerp * screenedLogTerm;
}

@nogc double customLAMPCoulombOmega23(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega2sFromLog(3, screenedLogTerm, rhoPerp);
}

@nogc double customLAMPCoulombOmega24(
    double screeningTemperature,
    double debyeTemperature,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double rhoPerp;
    double screenedLogTerm = customLAMPCoulombScreenedLogTerm(
        screeningTemperature, debyeTemperature, electronDensity, lampI, lampJ, rhoPerp);
    return customLAMPCoulombOmega2sFromLog(4, screenedLogTerm, rhoPerp);
}

@nogc bool customLAMPUsesEmpiricalNeutralPair(int i, int j)
{
    return (i == 0 && j == 4) || (i == 1 && j == 4) || i == 2 || (i == 3 && j == 4) || (i == 4 && j == 4);
}

@nogc bool customLAMPUsesEmpiricalIonNeutralPair(int neutralSpecies, int ionSpecies)
{
    return (neutralSpecies == 0 && ionSpecies == 7)
        || (neutralSpecies == 1 && ionSpecies == 8)
        || (neutralSpecies == 2 && ionSpecies == 5)
        || (neutralSpecies == 3 && ionSpecies == 6)
        || (neutralSpecies == 4 && ionSpecies == 9);
}

@nogc bool customLAMPUsesPhenomenologicalIonNeutralPair(int neutralSpecies, int ionSpecies)
{
    return (neutralSpecies == 0 && (ionSpecies == 5 || ionSpecies == 6))
        || (neutralSpecies == 1 && (ionSpecies == 5 || ionSpecies == 6))
        || (neutralSpecies == 2 && ionSpecies == 6)
        || (neutralSpecies == 3 && ionSpecies == 5);
}

@nogc double customLAMPCollisionOmega11(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega11(Te, i);
        return customLAMPCoulombOmega11(Te, Te, electronDensity, i, j);
    }
    if (i <= 4 && j <= 4) {
        if (customLAMPUsesEmpiricalNeutralPair(i, j)) {
            return customLAMPEmpiricalNeutralOmega11(Th, i, j);
        }
        return customLAMPPhenomenologicalNN(Th, i, j);
    }
    if (i <= 4 && j <= 9) {
        if (customLAMPUsesEmpiricalIonNeutralPair(i, j)) {
            return customLAMPEmpiricalIonNeutralOmega11(Th, i, j);
        }
        if (customLAMPUsesPhenomenologicalIonNeutralPair(i, j)) {
            return customLAMPPhenomenologicalNI(Th, i, j, false);
        }
        return customLAMPScaledPairOmega11(Th, i);
    }
    if (i >= 5 && j <= 9) return customLAMPCoulombOmega11(Th, Te, electronDensity, i, j);
    throw new Error("unsupported CustomLAMP collision pair");
}

@nogc double customLAMPCollisionOmega22(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega(Te, i, 22);
        return customLAMPCoulombOmega22(Te, Te, electronDensity, i, j);
    }
    if (i <= 4 && j <= 4) {
        if (customLAMPUsesEmpiricalNeutralPair(i, j)) {
            return customLAMPEmpiricalNeutralOmega22(Th, i, j);
        }
        return customLAMPPhenomenologicalNN22(Th, i, j);
    }
    if (i <= 4 && j <= 9) {
        if (customLAMPUsesEmpiricalIonNeutralPair(i, j)) {
            return customLAMPEmpiricalIonNeutralOmega22(Th, i, j);
        }
        if (customLAMPUsesPhenomenologicalIonNeutralPair(i, j)) {
            return customLAMPPhenomenologicalNI(Th, i, j, true);
        }
        return customLAMPScaledPairOmega22(Th, i);
    }
    if (i >= 5 && j <= 9) return customLAMPCoulombOmega22(Th, Te, electronDensity, i, j);
    throw new Error("unsupported CustomLAMP collision pair");
}

@nogc double customLAMPCollisionOmega12(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega(Te, i, 12);
        return customLAMPCoulombOmega12(Te, Te, electronDensity, i, j);
    }
    if (i <= 4 && j <= 4) {
        if (customLAMPUsesEmpiricalNeutralPair(i, j)) {
            return customLAMPEmpiricalNeutralOmega12(Th, i, j);
        }
        return customLAMPPhenomenologicalNN12(Th, i, j);
    }
    if (i <= 4 && j <= 9) {
        if (customLAMPUsesEmpiricalIonNeutralPair(i, j)) {
            return customLAMPEmpiricalIonNeutralOmega12(Th, i, j);
        }
        if (customLAMPUsesPhenomenologicalIonNeutralPair(i, j)) {
            return customLAMPPhenomenologicalNIByOrder(Th, i, j, 12);
        }
        return customLAMPScaledPairOmega12(Th, i);
    }
    if (i >= 5 && j <= 9) return customLAMPCoulombOmega12(Th, Te, electronDensity, i, j);
    throw new Error("unsupported CustomLAMP collision pair");
}

@nogc double customLAMPCollisionOmega13(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega(Te, i, 13);
        return customLAMPCoulombOmega13(Te, Te, electronDensity, i, j);
    }
    if (i <= 4 && j <= 4) {
        if (customLAMPUsesEmpiricalNeutralPair(i, j)) {
            return customLAMPEmpiricalNeutralOmega13(Th, i, j);
        }
        return customLAMPPhenomenologicalNN13(Th, i, j);
    }
    if (i <= 4 && j <= 9) {
        if (customLAMPUsesEmpiricalIonNeutralPair(i, j)) {
            return customLAMPEmpiricalIonNeutralOmega13(Th, i, j);
        }
        if (customLAMPUsesPhenomenologicalIonNeutralPair(i, j)) {
            return customLAMPPhenomenologicalNIByOrder(Th, i, j, 13);
        }
        return customLAMPScaledPairOmega13(Th, i);
    }
    if (i >= 5 && j <= 9) return customLAMPCoulombOmega13(Th, Te, electronDensity, i, j);
    throw new Error("unsupported CustomLAMP collision pair");
}

@nogc double customLAMPCollisionOmega14(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega(Te, i, 14);
        return customLAMPCoulombOmega14(Te, Te, electronDensity, i, j);
    }
    throw new Error("CustomLAMP Omega14 is only used for electron pairs");
}

@nogc double customLAMPCollisionOmega15(double Th, double Te, double electronDensity, int lampI, int lampJ)
{
    int i = lampI < lampJ ? lampI : lampJ;
    int j = lampI < lampJ ? lampJ : lampI;
    if (j == 10) {
        if (i <= 4) return customLAMPElectronNeutralOmega(Te, i, 15);
        return customLAMPCoulombOmega15(Te, Te, electronDensity, i, j);
    }
    throw new Error("CustomLAMP Omega15 is only used for electron pairs");
}

@nogc double customLAMPSpeciesTemperature(ref const(GasState) gs, int lampIndex)
{
    return lampIndex == 10 ? customLAMPReal(gs.T_modes[0]) : customLAMPReal(gs.T);
}

@nogc number customLAMPBinaryDiffusion(
    ref const(GasState) gs,
    double totalNumberDensity,
    double electronDensity,
    int lampI,
    int lampJ)
{
    double massI = customLAMPSpeciesMass(lampI);
    double massJ = customLAMPSpeciesMass(lampJ);
    double Ti = customLAMPSpeciesTemperature(gs, lampI);
    double Tj = customLAMPSpeciesTemperature(gs, lampJ);
    double Th = customLAMPReal(gs.T);
    double Te = customLAMPReal(gs.T_modes[0]);
    double omega11 = customLAMPCollisionOmega11(Th, Te, electronDensity, lampI, lampJ);
    if (omega11 <= 0.0) throw new Error("CustomLAMP collision integral must be positive");
    double reducedMass = massI * massJ / (massI + massJ);
    double effectiveTemperature = 1.0 / ((massI / Tj + massJ / Ti) / (massI + massJ));
    double thetaJI = Tj / Ti;
    double gIJ = pow(massJ + massI * thetaJI, 2.0)
        * pow(massI + massJ * thetaJI, 3.0)
        / (pow(thetaJI, 2.5) * pow(massI + massJ, 5.0));
    return customLAMPNumberFromDouble(
        3.0 * Ti * gIJ / (16.0 * totalNumberDensity * omega11)
        * sqrt(2.0 * CUSTOM_LAMP_PI * CUSTOM_LAMP_KB / (reducedMass * effectiveTemperature))
    );
}

@nogc void customLAMPCopyHeavyNumberDensities(ref const(number[]) lampNumberDensities, ref number[] heavyNumberDensities)
{
    debug { assert(heavyNumberDensities.length == CUSTOM_LAMP_HEAVY_SPECIES); }
    foreach (index; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        heavyNumberDensities[index] = lampNumberDensities[index];
    }
}

@nogc void customLAMPBuildHeavyCollisionMatrices(
    ref const(GasState) gs,
    double electronDensity,
    ref number[][] collision11,
    ref number[][] collision12,
    ref number[][] collision13,
    ref number[][] collision22)
{
    debug {
        assert(collision11.length == CUSTOM_LAMP_HEAVY_SPECIES);
        assert(collision12.length == CUSTOM_LAMP_HEAVY_SPECIES);
        assert(collision13.length == CUSTOM_LAMP_HEAVY_SPECIES);
        assert(collision22.length == CUSTOM_LAMP_HEAVY_SPECIES);
    }
    double Th = customLAMPReal(gs.T);
    double Te = customLAMPReal(gs.T_modes[0]);
    foreach (i; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        foreach (j; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            collision11[i][j] = customLAMPCollisionOmega11(Th, Te, electronDensity, to!int(i), to!int(j));
            collision12[i][j] = customLAMPCollisionOmega12(Th, Te, electronDensity, to!int(i), to!int(j));
            collision13[i][j] = customLAMPCollisionOmega13(Th, Te, electronDensity, to!int(i), to!int(j));
            collision22[i][j] = customLAMPCollisionOmega22(Th, Te, electronDensity, to!int(i), to!int(j));
        }
    }
}

@nogc void customLAMPBuildViscositySystem(
    double heavyTemperature,
    ref const(number[]) zfHeavy,
    ref const(number[][]) collision11,
    ref const(number[][]) collision22,
    ref number[][] core,
    ref number[] borderCol,
    ref number[] borderRow)
{
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        double massRow = customLAMPSpeciesMass(to!int(row));
        foreach (col; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            double massCol = customLAMPSpeciesMass(to!int(col));
            number value = 0.0;
            foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
                double massSpecies = customLAMPSpeciesMass(to!int(species));
                double deltaRowCol = row == col ? 1.0 : 0.0;
                double deltaColSpecies = col == species ? 1.0 : 0.0;
                value += zfHeavy[species] * sqrt(massSpecies) / pow(massRow + massSpecies, 1.5)
                    * (
                        10.0 / 3.0
                        * collision11[row][species]
                        * (deltaRowCol - deltaColSpecies)
                        * massCol
                        + 2.0
                        * massSpecies
                        * collision22[row][species]
                        * (deltaRowCol + deltaColSpecies)
                    );
            }
            core[row][col] = value * (8.0 * zfHeavy[row] * massRow / massCol);
        }
        borderCol[row] = zfHeavy[row] * sqrt(massRow);
        borderRow[row] = zfHeavy[row];
    }
}

@nogc void customLAMPLUFactorInPlace(ref number[][] values, ref int[] pivots)
{
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) pivots[row] = to!int(row);
    foreach (col; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        size_t pivotRow = col;
        number pivotAbs = abs(values[col][col]);
        foreach (row; col + 1 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            number candidateAbs = abs(values[row][col]);
            if (candidateAbs > pivotAbs) {
                pivotAbs = candidateAbs;
                pivotRow = row;
            }
        }
        if (pivotRow != col) {
            foreach (k; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
                number tmp = values[col][k];
                values[col][k] = values[pivotRow][k];
                values[pivotRow][k] = tmp;
            }
            int tmpPivot = pivots[col];
            pivots[col] = pivots[pivotRow];
            pivots[pivotRow] = tmpPivot;
        }
        number pivot = values[col][col];
        if (abs(pivot) <= CUSTOM_LAMP_PIVOT_FLOOR) continue;
        foreach (row; col + 1 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            values[row][col] /= pivot;
            number multiplier = values[row][col];
            if (multiplier == 0.0) continue;
            foreach (k; col + 1 .. CUSTOM_LAMP_HEAVY_SPECIES) {
                values[row][k] -= multiplier * values[col][k];
            }
        }
    }
}

@nogc void customLAMPSolveLU(
    ref number[][] lu,
    ref const(int[]) pivots,
    ref const(number[]) rhs,
    ref number[] scratch,
    ref number[] solution)
{
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        scratch[row] = rhs[pivots[row]];
    }
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        foreach (col; 0 .. row) {
            scratch[row] -= lu[row][col] * scratch[col];
        }
    }
    foreach_reverse (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        number acc = scratch[row];
        foreach (col; row + 1 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            acc -= lu[row][col] * solution[col];
        }
        number pivot = lu[row][row];
        if (abs(pivot) <= CUSTOM_LAMP_PIVOT_FLOOR) {
            throw new Error("singular matrix in CustomLAMP viscosity solve");
        }
        solution[row] = acc / pivot;
    }
}

@nogc number customLAMPBorderedDeterminantRatio(
    ref const(number[][]) core,
    ref const(number[]) borderCol,
    ref const(number[]) borderRow,
    ref number[][] lu,
    ref int[] pivots,
    ref number[] scratch,
    ref number[] solution)
{
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        foreach (col; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            lu[row][col] = core[row][col];
        }
    }
    customLAMPLUFactorInPlace(lu, pivots);
    customLAMPSolveLU(lu, pivots, borderCol, scratch, solution);
    number dot = 0.0;
    foreach (index; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        dot += borderRow[index] * solution[index];
    }
    return -dot;
}

@nogc number customLAMPViscosity(
    ref const(GasState) gs,
    ref const(number[]) heavyNumberDensities,
    ref const(number[][]) collision11,
    ref const(number[][]) collision22,
    ref number[][] core,
    ref number[] borderCol,
    ref number[] borderRow,
    ref number[][] lu,
    ref int[] pivots,
    ref number[] scratch,
    ref number[] solution)
{
    customLAMPBuildViscositySystem(customLAMPReal(gs.T), heavyNumberDensities, collision11, collision22, core, borderCol, borderRow);
    number ratio = customLAMPBorderedDeterminantRatio(core, borderCol, borderRow, lu, pivots, scratch, solution);
    return -2.5 * sqrt(2.0 * CUSTOM_LAMP_PI * CUSTOM_LAMP_KB) * sqrt(gs.T) * ratio;
}

@nogc void customLAMPLUFactorInPlaceN(ref number[][] values, ref int[] pivots, int size)
{
    foreach (row; 0 .. size) pivots[row] = row;
    foreach (col; 0 .. size) {
        int pivotRow = col;
        number pivotAbs = abs(values[col][col]);
        foreach (row; col + 1 .. size) {
            number candidateAbs = abs(values[row][col]);
            if (candidateAbs > pivotAbs) {
                pivotAbs = candidateAbs;
                pivotRow = row;
            }
        }
        if (pivotRow != col) {
            foreach (k; 0 .. size) {
                number tmp = values[col][k];
                values[col][k] = values[pivotRow][k];
                values[pivotRow][k] = tmp;
            }
            int tmpPivot = pivots[col];
            pivots[col] = pivots[pivotRow];
            pivots[pivotRow] = tmpPivot;
        }
        number pivot = values[col][col];
        if (abs(pivot) <= CUSTOM_LAMP_PIVOT_FLOOR) continue;
        foreach (row; col + 1 .. size) {
            values[row][col] /= pivot;
            number multiplier = values[row][col];
            if (multiplier == 0.0) continue;
            foreach (k; col + 1 .. size) {
                values[row][k] -= multiplier * values[col][k];
            }
        }
    }
}

@nogc void customLAMPSolveLUN(
    ref number[][] lu,
    ref const(int[]) pivots,
    ref const(number[]) rhs,
    ref number[] scratch,
    ref number[] solution,
    int size)
{
    foreach (row; 0 .. size) {
        scratch[row] = rhs[pivots[row]];
        solution[row] = 0.0;
    }
    foreach (row; 0 .. size) {
        foreach (col; 0 .. row) {
            scratch[row] -= lu[row][col] * scratch[col];
        }
    }
    foreach_reverse (row; 0 .. size) {
        number acc = scratch[row];
        foreach (col; row + 1 .. size) {
            acc -= lu[row][col] * solution[col];
        }
        number pivot = lu[row][row];
        if (abs(pivot) <= CUSTOM_LAMP_PIVOT_FLOOR) {
            throw new Error("singular matrix in CustomLAMP thermal-conductivity solve");
        }
        solution[row] = acc / pivot;
    }
}

@nogc number customLAMPBorderedDeterminantRatioN(
    ref const(number[][]) core,
    ref const(number[]) borderCol,
    ref const(number[]) borderRow,
    ref number[][] lu,
    ref int[] pivots,
    ref number[] scratch,
    ref number[] solution,
    int size)
{
    foreach (row; 0 .. size) {
        foreach (col; 0 .. size) {
            lu[row][col] = core[row][col];
        }
    }
    customLAMPLUFactorInPlaceN(lu, pivots, size);
    customLAMPSolveLUN(lu, pivots, borderCol, scratch, solution, size);
    number dot = 0.0;
    foreach (index; 0 .. size) {
        dot += borderRow[index] * solution[index];
    }
    return -dot;
}

@nogc void customLAMPBuildHeavyThermalSystem(
    ref const(number[]) zfHeavy,
    ref const(number[][]) collision11,
    ref const(number[][]) collision12,
    ref const(number[][]) collision13,
    ref const(number[][]) collision22,
    ref number[][] core,
    ref number[] borderCol,
    ref number[] borderRow)
{
    foreach (index; 0 .. CUSTOM_LAMP_THERMAL_SYSTEM) {
        borderCol[index] = 0.0;
        borderRow[index] = 0.0;
    }
    foreach (row; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        double massRow = customLAMPSpeciesMass(to!int(row));
        foreach (col; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            double massCol = customLAMPSpeciesMass(to!int(col));
            number q00 = 0.0;
            number q01 = 0.0;
            number q11 = 0.0;
            foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
                double massSpecies = customLAMPSpeciesMass(to!int(species));
                double deltaRowCol = row == col ? 1.0 : 0.0;
                double deltaColSpecies = col == species ? 1.0 : 0.0;
                double deltaRowSpecies = row == species ? 1.0 : 0.0;
                q00 += zfHeavy[species] * sqrt(massRow) / sqrt(massRow + massSpecies)
                    * collision11[row][species]
                    * (
                        zfHeavy[row] * sqrt(massSpecies) / sqrt(massCol)
                        * (deltaRowCol - deltaColSpecies)
                        - zfHeavy[col] * sqrt(massSpecies * massCol) / massRow
                        * (1.0 - deltaRowSpecies)
                    );
                q01 += zfHeavy[species] * pow(massSpecies / (massRow + massSpecies), 1.5)
                    * (2.5 * collision11[row][species] - 3.0 * collision12[row][species])
                    * (deltaRowCol - deltaColSpecies);
                q11 += zfHeavy[species] * sqrt(massSpecies) / pow(massRow + massSpecies, 2.5)
                    * (
                        (deltaRowCol - deltaColSpecies)
                        * (
                            1.25 * (
                                6.0 * massCol * massCol
                                + 5.0 * massSpecies * massSpecies
                            ) * collision11[row][species]
                            - 15.0 * massSpecies * massSpecies * collision12[row][species]
                            + 12.0 * massSpecies * massSpecies * collision13[row][species]
                        )
                        + (deltaRowCol + deltaColSpecies)
                        * 4.0 * massCol * massSpecies * collision22[row][species]
                    );
            }
            q00 *= 8.0;
            q01 *= 8.0 * zfHeavy[row] * pow(massRow / massCol, 1.5);
            q11 *= 8.0 * zfHeavy[row] * pow(massRow / massCol, 1.5);
            number q10 = massCol / massRow * q01;
            core[row][col] = q00;
            core[row][col + CUSTOM_LAMP_HEAVY_SPECIES] = q01;
            core[row + CUSTOM_LAMP_HEAVY_SPECIES][col] = q10;
            core[row + CUSTOM_LAMP_HEAVY_SPECIES][col + CUSTOM_LAMP_HEAVY_SPECIES] = q11;
        }
    }
    foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        double massSpecies = customLAMPSpeciesMass(to!int(species));
        borderCol[CUSTOM_LAMP_HEAVY_SPECIES + species] = zfHeavy[species];
        borderRow[CUSTOM_LAMP_HEAVY_SPECIES + species] = zfHeavy[species] / sqrt(massSpecies);
    }
}

@nogc number customLAMPHeavyTranslationalThermalConductivity(
    ref const(GasState) gs,
    ref const(number[]) heavyNumberDensities,
    ref const(number[][]) collision11,
    ref const(number[][]) collision12,
    ref const(number[][]) collision13,
    ref const(number[][]) collision22,
    ref number[][] core,
    ref number[] borderCol,
    ref number[] borderRow,
    ref number[][] lu,
    ref int[] pivots,
    ref number[] scratch,
    ref number[] solution)
{
    customLAMPBuildHeavyThermalSystem(
        heavyNumberDensities,
        collision11,
        collision12,
        collision13,
        collision22,
        core,
        borderCol,
        borderRow
    );
    number ratio = customLAMPBorderedDeterminantRatioN(
        core,
        borderCol,
        borderRow,
        lu,
        pivots,
        scratch,
        solution,
        CUSTOM_LAMP_THERMAL_SYSTEM
    );
    number cht = -75.0 * CUSTOM_LAMP_KB * sqrt(2.0 * CUSTOM_LAMP_PI * CUSTOM_LAMP_KB) / 8.0;
    return cht * sqrt(gs.T) * ratio;
}

@nogc void customLAMPBuildHeavyBinaryDiffusion(
    ref const(GasState) gs,
    double totalNumberDensity,
    double electronDensity,
    ref number[][] binaryDiffusion)
{
    foreach (i; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        foreach (j; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            binaryDiffusion[i][j] = customLAMPBinaryDiffusion(
                gs,
                totalNumberDensity,
                electronDensity,
                to!int(i),
                to!int(j)
            );
        }
    }
}

@nogc double customLAMPRotationalHeatCapacity(int lampIndex)
{
    if (lampIndex == 0 || lampIndex == 1 || lampIndex == 4 || lampIndex == 7 || lampIndex == 8 || lampIndex == 9) {
        return CUSTOM_LAMP_KB / customLAMPSpeciesMass(lampIndex);
    }
    return 0.0;
}

@nogc number customLAMPHeavyInternalRotationalThermalConductivity(
    ref const(GasState) gs,
    ref const(number[]) lampNumberDensities,
    double totalNumberDensity,
    ref const(number[][]) binaryDiffusion)
{
    number conductivity = 0.0;
    foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        number diffusionDenominator = 0.0;
        foreach (partner; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            number partnerD = fmax(binaryDiffusion[species][partner], CUSTOM_LAMP_INTERNAL_DENOM_FLOOR);
            diffusionDenominator += (lampNumberDensities[partner] / totalNumberDensity) / partnerD;
        }
        diffusionDenominator = fmax(diffusionDenominator, CUSTOM_LAMP_INTERNAL_DENOM_FLOOR);
        conductivity += gs.p / (CUSTOM_LAMP_KB * gs.T)
            * (lampNumberDensities[species] / totalNumberDensity)
            * customLAMPRotationalHeatCapacity(to!int(species))
            * customLAMPSpeciesMass(to!int(species))
            / diffusionDenominator;
    }
    return conductivity;
}

@nogc bool customLAMPIsMolecular(int lampIndex)
{
    return lampIndex == 0 || lampIndex == 1 || lampIndex == 4
        || lampIndex == 7 || lampIndex == 8 || lampIndex == 9;
}

@nogc double customLAMPVibrationalModeEnthalpy(
    double electronTemperature,
    double waveNumberCm,
    double anharmonicityCm,
    double cutoffEnergyEv)
{
    double partitionSum = 0.0;
    double energySum = 0.0;
    foreach (level; 0 .. 101) {
        double x = to!double(level) + 0.5;
        double vibrationalEnergyEv = (waveNumberCm * x - anharmonicityCm * x * x)
            * CUSTOM_LAMP_H / CUSTOM_LAMP_EV2J * (100.0 * CUSTOM_LAMP_C0);
        if (vibrationalEnergyEv < cutoffEnergyEv) {
            double boltzmannFactor = exp(-vibrationalEnergyEv * CUSTOM_LAMP_EV2J
                / (CUSTOM_LAMP_KB * electronTemperature));
            partitionSum += boltzmannFactor;
            energySum += vibrationalEnergyEv * CUSTOM_LAMP_EV2J * boltzmannFactor;
        }
    }
    return energySum / partitionSum;
}

@nogc double customLAMPMolecularVibrationalEnthalpy(double heavyTemperature, int lampIndex)
{
    final switch (lampIndex) {
    case 0: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 2358.57, 14.324, 9.759);
    case 1: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 1580.193, 11.981, 5.115);
    case 2: return 0.0;
    case 3: return 0.0;
    case 4: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 1904.20, 14.075, 6.497);
    case 5: return 0.0;
    case 6: return 0.0;
    case 7: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 2207.00, 16.10, 8.713);
    case 8: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 1904.77, 16.259, 6.663);
    case 9: return customLAMPVibrationalModeEnthalpy(heavyTemperature, 2376.42, 16.262, 10.851);
    case 10: return 0.0;
    }
}

@nogc double customLAMPAtomicElectronicLevelEnergy(int lampIndex, int level)
{
    if (lampIndex == 2) {
        final switch (level) {
        case 0: return 0.0 * CUSTOM_LAMP_EC;
        case 1: return 2.384 * CUSTOM_LAMP_EC;
        case 2: return 3.566 * CUSTOM_LAMP_EC;
        case 3: return 10.331 * CUSTOM_LAMP_EC;
        case 4: return 10.685 * CUSTOM_LAMP_EC;
        case 5: return 10.929 * CUSTOM_LAMP_EC;
        case 6: return 11.603 * CUSTOM_LAMP_EC;
        case 7: return 11.756 * CUSTOM_LAMP_EC;
        case 8: return 11.841 * CUSTOM_LAMP_EC;
        case 9: return 11.996 * CUSTOM_LAMP_EC;
        case 10: return 12.005 * CUSTOM_LAMP_EC;
        case 11: return 12.124 * CUSTOM_LAMP_EC;
        case 12: return 12.356 * CUSTOM_LAMP_EC;
        case 13: return 12.854 * CUSTOM_LAMP_EC;
        case 14: return 12.917 * CUSTOM_LAMP_EC;
        }
    } else if (lampIndex == 3) {
        final switch (level) {
        case 0: return 0.0159 * CUSTOM_LAMP_EC;
        case 1: return 1.967 * CUSTOM_LAMP_EC;
        case 2: return 4.190 * CUSTOM_LAMP_EC;
        case 3: return 9.146 * CUSTOM_LAMP_EC;
        case 4: return 9.521 * CUSTOM_LAMP_EC;
        case 5: return 10.740 * CUSTOM_LAMP_EC;
        case 6: return 10.989 * CUSTOM_LAMP_EC;
        case 7: return 11.838 * CUSTOM_LAMP_EC;
        case 8: return 11.930 * CUSTOM_LAMP_EC;
        case 9: return 12.079 * CUSTOM_LAMP_EC;
        case 10: return 12.087 * CUSTOM_LAMP_EC;
        case 11: return 12.286 * CUSTOM_LAMP_EC;
        case 12: return 12.359 * CUSTOM_LAMP_EC;
        case 13: return 12.541 * CUSTOM_LAMP_EC;
        case 14: return 12.661 * CUSTOM_LAMP_EC;
        }
    } else if (lampIndex == 5) {
        final switch (level) {
        case 0: return 0.007 * CUSTOM_LAMP_EC;
        case 1: return 1.899 * CUSTOM_LAMP_EC;
        case 2: return 4.053 * CUSTOM_LAMP_EC;
        case 3: return 5.724 * CUSTOM_LAMP_EC;
        case 4: return 11.437 * CUSTOM_LAMP_EC;
        case 5: return 13.542 * CUSTOM_LAMP_EC;
        case 6: return 17.471 * CUSTOM_LAMP_EC;
        case 7: return 18.471 * CUSTOM_LAMP_EC;
        case 8: return 18.497 * CUSTOM_LAMP_EC;
        case 9: return 19.234 * CUSTOM_LAMP_EC;
        case 10: return 20.409 * CUSTOM_LAMP_EC;
        case 11: return 20.655 * CUSTOM_LAMP_EC;
        case 12: return 20.676 * CUSTOM_LAMP_EC;
        case 13: return 20.940 * CUSTOM_LAMP_EC;
        case 14: return 21.154 * CUSTOM_LAMP_EC;
        }
    } else if (lampIndex == 6) {
        final switch (level) {
        case 0: return 0.0 * CUSTOM_LAMP_EC;
        case 1: return 3.325 * CUSTOM_LAMP_EC;
        case 2: return 5.018 * CUSTOM_LAMP_EC;
        case 3: return 14.875 * CUSTOM_LAMP_EC;
        case 4: return 20.580 * CUSTOM_LAMP_EC;
        case 5: return 22.982 * CUSTOM_LAMP_EC;
        case 6: return 23.430 * CUSTOM_LAMP_EC;
        case 7: return 24.265 * CUSTOM_LAMP_EC;
        case 8: return 25.286 * CUSTOM_LAMP_EC;
        case 9: return 25.646 * CUSTOM_LAMP_EC;
        case 10: return 25.661 * CUSTOM_LAMP_EC;
        case 11: return 25.839 * CUSTOM_LAMP_EC;
        case 12: return 26.237 * CUSTOM_LAMP_EC;
        case 13: return 26.305 * CUSTOM_LAMP_EC;
        case 14: return 26.369 * CUSTOM_LAMP_EC;
        }
    }
    return 0.0;
}

@nogc double customLAMPAtomicElectronicLevelDegeneracy(int lampIndex, int level)
{
    if (lampIndex == 2) {
        final switch (level) {
        case 0: return 4.0;
        case 1: return 10.0;
        case 2: return 6.0;
        case 3: return 12.0;
        case 4: return 6.0;
        case 5: return 12.0;
        case 6: return 2.0;
        case 7: return 20.0;
        case 8: return 12.0;
        case 9: return 4.0;
        case 10: return 10.0;
        case 11: return 6.0;
        case 12: return 10.0;
        case 13: return 12.0;
        case 14: return 6.0;
        }
    } else if (lampIndex == 3) {
        final switch (level) {
        case 0: return 9.0;
        case 1: return 5.0;
        case 2: return 1.0;
        case 3: return 5.0;
        case 4: return 3.0;
        case 5: return 15.0;
        case 6: return 9.0;
        case 7: return 5.0;
        case 8: return 3.0;
        case 9: return 25.0;
        case 10: return 15.0;
        case 11: return 15.0;
        case 12: return 9.0;
        case 13: return 15.0;
        case 14: return 5.0;
        }
    } else if (lampIndex == 5) {
        final switch (level) {
        case 0: return 9.0;
        case 1: return 5.0;
        case 2: return 1.0;
        case 3: return 5.0;
        case 4: return 15.0;
        case 5: return 9.0;
        case 6: return 5.0;
        case 7: return 9.0;
        case 8: return 3.0;
        case 9: return 3.0;
        case 10: return 3.0;
        case 11: return 15.0;
        case 12: return 3.0;
        case 13: return 3.0;
        case 14: return 9.0;
        }
    } else if (lampIndex == 6) {
        final switch (level) {
        case 0: return 4.0;
        case 1: return 10.0;
        case 2: return 6.0;
        case 3: return 12.0;
        case 4: return 10.0;
        case 5: return 12.0;
        case 6: return 6.0;
        case 7: return 2.0;
        case 8: return 2.0;
        case 9: return 20.0;
        case 10: return 10.0;
        case 11: return 12.0;
        case 12: return 10.0;
        case 13: return 4.0;
        case 14: return 6.0;
        }
    }
    return 0.0;
}

@nogc double customLAMPAtomicElectronicEnthalpy(double electronTemperature, int lampIndex)
{
    double partitionSum = 0.0;
    double energySum = 0.0;
    foreach (level; 0 .. 15) {
        double energy = customLAMPAtomicElectronicLevelEnergy(lampIndex, to!int(level));
        double degeneracy = customLAMPAtomicElectronicLevelDegeneracy(lampIndex, to!int(level));
        double boltzmannFactor = exp(-energy / (CUSTOM_LAMP_KB * electronTemperature));
        partitionSum += degeneracy * boltzmannFactor;
        energySum += energy * degeneracy * boltzmannFactor;
    }
    return energySum / partitionSum;
}

@nogc void customLAMPBuildHeavyModalHeatCapacities(
    ref const(GasState) gs,
    ref number[] vibrationalHeatCapacity,
    ref number[] electronicHeatCapacity)
{
    debug {
        assert(vibrationalHeatCapacity.length == CUSTOM_LAMP_HEAVY_SPECIES);
        assert(electronicHeatCapacity.length == CUSTOM_LAMP_HEAVY_SPECIES);
    }
    double heavyTemperature = customLAMPReal(gs.T);
    double electronTemperature = customLAMPReal(gs.T_modes[0]);
    foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        int lampIndex = to!int(species);
        vibrationalHeatCapacity[species] = 0.0;
        electronicHeatCapacity[species] = 0.0;
        double mass = customLAMPSpeciesMass(lampIndex);
        if (customLAMPIsMolecular(lampIndex)) {
            double high = customLAMPMolecularVibrationalEnthalpy(heavyTemperature + 100.0, lampIndex) / mass;
            double base = customLAMPMolecularVibrationalEnthalpy(heavyTemperature, lampIndex) / mass;
            double modalHeatCapacity = fmax((high - base) / 100.0, 0.0);
            vibrationalHeatCapacity[species] = fmax(modalHeatCapacity - customLAMPRotationalHeatCapacity(lampIndex), 0.0);
        } else if (lampIndex == 2 || lampIndex == 3 || lampIndex == 5 || lampIndex == 6) {
            double high = customLAMPAtomicElectronicEnthalpy(electronTemperature + 100.0, lampIndex) / mass;
            double base = customLAMPAtomicElectronicEnthalpy(electronTemperature, lampIndex) / mass;
            electronicHeatCapacity[species] = fmax((high - base) / 100.0, 0.0);
        }
    }
}

@nogc number customLAMPHeavyInternalThermalConductivityFromModalHeatCapacity(
    ref const(GasState) gs,
    ref const(number[]) lampNumberDensities,
    double totalNumberDensity,
    ref const(number[]) modalHeatCapacity,
    ref const(number[][]) binaryDiffusion)
{
    number conductivity = 0.0;
    foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        number diffusionDenominator = 0.0;
        foreach (partner; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
            number partnerD = fmax(binaryDiffusion[species][partner], CUSTOM_LAMP_INTERNAL_DENOM_FLOOR);
            diffusionDenominator += (lampNumberDensities[partner] / totalNumberDensity) / partnerD;
        }
        diffusionDenominator = fmax(diffusionDenominator, CUSTOM_LAMP_INTERNAL_DENOM_FLOOR);
        conductivity += gs.p / (CUSTOM_LAMP_KB * gs.T)
            * (lampNumberDensities[species] / totalNumberDensity)
            * modalHeatCapacity[species]
            * customLAMPSpeciesMass(to!int(species))
            / diffusionDenominator;
    }
    return conductivity;
}

@nogc void customLAMPDevotoElectronQ(
    ref const(GasState) gs,
    ref const(number[]) lampNumberDensities,
    double electronDensity,
    out number q11,
    out number q12,
    out number q22)
{
    q11 = 0.0;
    q12 = 0.0;
    q22 = 0.0;
    double Th = customLAMPReal(gs.T);
    double Te = customLAMPReal(gs.T_modes[0]);
    foreach (species; 0 .. CUSTOM_LAMP_HEAVY_SPECIES) {
        double omega11 = customLAMPCollisionOmega11(Th, Te, electronDensity, 10, to!int(species));
        double omega12 = customLAMPCollisionOmega12(Th, Te, electronDensity, 10, to!int(species));
        double omega13 = customLAMPCollisionOmega13(Th, Te, electronDensity, 10, to!int(species));
        double omega14 = customLAMPCollisionOmega14(Th, Te, electronDensity, 10, to!int(species));
        double omega15 = customLAMPCollisionOmega15(Th, Te, electronDensity, 10, to!int(species));
        q11 += lampNumberDensities[species] * (
            25.0 / 4.0 * omega11
            - 15.0 * omega12
            + 12.0 * omega13
        );
        q12 += lampNumberDensities[species] * (
            175.0 / 16.0 * omega11
            - 315.0 / 8.0 * omega12
            + 57.0 * omega13
            - 30.0 * omega14
        );
        q22 += lampNumberDensities[species] * (
            1225.0 / 64.0 * omega11
            - 735.0 / 8.0 * omega12
            + 399.0 / 2.0 * omega13
            - 210.0 * omega14
            + 90.0 * omega15
        );
    }
    q11 *= 8.0 * electronDensity;
    q12 *= 8.0 * electronDensity;
    q22 *= 8.0 * electronDensity;

    double electronOmega22 = customLAMPCoulombOmega22(Te, Te, electronDensity, 10, 10);
    double electronOmega23 = customLAMPCoulombOmega23(Te, Te, electronDensity, 10, 10);
    double electronOmega24 = customLAMPCoulombOmega24(Te, Te, electronDensity, 10, 10);
    double eePrefactor = 8.0 * sqrt(2.0) * electronDensity * electronDensity;
    q11 += eePrefactor * electronOmega22;
    q12 += eePrefactor * (7.0 / 4.0 * electronOmega22 - 2.0 * electronOmega23);
    q22 += eePrefactor * (77.0 / 16.0 * electronOmega22 - 7.0 * electronOmega23 + 5.0 * electronOmega24);
}

@nogc number customLAMPElectronTranslationalThermalConductivity(
    ref const(GasState) gs,
    ref const(number[]) lampNumberDensities,
    double electronDensity)
{
    number q11, q12, q22;
    customLAMPDevotoElectronQ(gs, lampNumberDensities, electronDensity, q11, q12, q22);
    number denominator = q11 - q12 * q12 / q22;
    number cet = 75.0 * CUSTOM_LAMP_KB * sqrt(2.0 * CUSTOM_LAMP_PI * CUSTOM_LAMP_KB / CUSTOM_LAMP_ME) / 8.0;
    number conductivity = cet * electronDensity * electronDensity * sqrt(gs.T_modes[0]) / denominator;
    return customLAMPReal(conductivity) < 0.0 ? customLAMPNumberFromDouble(0.0) : conductivity;
}

class CustomLAMPTransProps : TransportPropertiesModel {
public:
    this(lua_State *L, string[] speciesNames)
    {
        mNSpecies = to!int(speciesNames.length);
        mGDTKToLAMPIndex.length = mNSpecies;
        mLampNumberDensities.length = CUSTOM_LAMP_N_SPECIES;
        mHeavyNumberDensities.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mVibrationalHeatCapacity.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mElectronicHeatCapacity.length = CUSTOM_LAMP_HEAVY_SPECIES;
        allocateMatrix(mCollision11, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mCollision12, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mCollision13, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mCollision22, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mHeavyBinaryDiffusion, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mViscosityCore, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mViscosityLU, CUSTOM_LAMP_HEAVY_SPECIES);
        allocateMatrix(mThermalCore, CUSTOM_LAMP_THERMAL_SYSTEM);
        allocateMatrix(mThermalLU, CUSTOM_LAMP_THERMAL_SYSTEM);
        mViscosityBorderCol.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mViscosityBorderRow.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mViscosityScratch.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mViscositySolution.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mViscosityPivots.length = CUSTOM_LAMP_HEAVY_SPECIES;
        mThermalBorderCol.length = CUSTOM_LAMP_THERMAL_SYSTEM;
        mThermalBorderRow.length = CUSTOM_LAMP_THERMAL_SYSTEM;
        mThermalScratch.length = CUSTOM_LAMP_THERMAL_SYSTEM;
        mThermalSolution.length = CUSTOM_LAMP_THERMAL_SYSTEM;
        mThermalPivots.length = CUSTOM_LAMP_THERMAL_SYSTEM;
        foreach (isp, speciesName; speciesNames) {
            int lampIndex = customLAMPLampSpeciesIndex(speciesName);
            if (lampIndex < 0) {
                throw new Error(
                    format("Species '%s' is not supported by CustomLAMP transport_model.", speciesName)
                );
            }
            mGDTKToLAMPIndex[isp] = lampIndex;
        }
    }

    @nogc
    void numberDensities(ref const(GasState) gs, ref number[] result)
    {
        debug { assert(result.length == CUSTOM_LAMP_N_SPECIES); }
        result[] = customLAMPNumberFromDouble(1.0e-10);
        foreach (gdtkIndex; 0 .. mNSpecies) {
            int lampIndex = mGDTKToLAMPIndex[gdtkIndex];
            result[lampIndex] = customLAMPNumberDensity(gs, gdtkIndex, lampIndex);
        }
    }

    @nogc
    override void updateTransProps(ref GasState gs, GasModel gm)
    {
        numberDensities(gs, mLampNumberDensities);
        double totalNumberDensity = 0.0;
        foreach (numberDensity; mLampNumberDensities) totalNumberDensity += customLAMPReal(numberDensity);
        double electronDensity = customLAMPReal(mLampNumberDensities[10]);
        customLAMPCopyHeavyNumberDensities(mLampNumberDensities, mHeavyNumberDensities);
        customLAMPBuildHeavyCollisionMatrices(gs, electronDensity, mCollision11, mCollision12, mCollision13, mCollision22);
        gs.mu = customLAMPViscosity(
            gs,
            mHeavyNumberDensities,
            mCollision11,
            mCollision22,
            mViscosityCore,
            mViscosityBorderCol,
            mViscosityBorderRow,
            mViscosityLU,
            mViscosityPivots,
            mViscosityScratch,
            mViscositySolution
        );
        customLAMPBuildHeavyBinaryDiffusion(gs, totalNumberDensity, electronDensity, mHeavyBinaryDiffusion);
        gs.k = customLAMPHeavyTranslationalThermalConductivity(
            gs,
            mHeavyNumberDensities,
            mCollision11,
            mCollision12,
            mCollision13,
            mCollision22,
            mThermalCore,
            mThermalBorderCol,
            mThermalBorderRow,
            mThermalLU,
            mThermalPivots,
            mThermalScratch,
            mThermalSolution
        ) + customLAMPHeavyInternalRotationalThermalConductivity(
            gs,
            mLampNumberDensities,
            totalNumberDensity,
            mHeavyBinaryDiffusion
        );
        customLAMPBuildHeavyModalHeatCapacities(gs, mVibrationalHeatCapacity, mElectronicHeatCapacity);
        if (gs.k_modes.length > 0) {
            gs.k_modes[0] = customLAMPHeavyInternalThermalConductivityFromModalHeatCapacity(
                gs,
                mLampNumberDensities,
                totalNumberDensity,
                mVibrationalHeatCapacity,
                mHeavyBinaryDiffusion
            ) + customLAMPHeavyInternalThermalConductivityFromModalHeatCapacity(
                gs,
                mLampNumberDensities,
                totalNumberDensity,
                mElectronicHeatCapacity,
                mHeavyBinaryDiffusion
            ) + customLAMPElectronTranslationalThermalConductivity(
                gs,
                mLampNumberDensities,
                electronDensity
            );
        }
    }

    @nogc
    override void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D)
    {
        numberDensities(gs, mLampNumberDensities);
        double totalNumberDensity = 0.0;
        foreach (numberDensity; mLampNumberDensities) totalNumberDensity += customLAMPReal(numberDensity);
        double electronDensity = customLAMPReal(mLampNumberDensities[10]);
        foreach (isp; 0 .. mNSpecies) {
            foreach (jsp; 0 .. mNSpecies) {
                if (isp == jsp) {
                    D[isp][jsp] = 0.0;
                } else {
                    D[isp][jsp] = customLAMPBinaryDiffusion(
                        gs,
                        totalNumberDensity,
                        electronDensity,
                        mGDTKToLAMPIndex[isp],
                        mGDTKToLAMPIndex[jsp]
                    );
                }
            }
        }
    }

private:
    int mNSpecies;
    int[] mGDTKToLAMPIndex;
    number[] mLampNumberDensities;
    number[] mHeavyNumberDensities;
    number[] mVibrationalHeatCapacity;
    number[] mElectronicHeatCapacity;
    number[][] mCollision11;
    number[][] mCollision12;
    number[][] mCollision13;
    number[][] mCollision22;
    number[][] mHeavyBinaryDiffusion;
    number[][] mViscosityCore;
    number[][] mViscosityLU;
    number[][] mThermalCore;
    number[][] mThermalLU;
    number[] mViscosityBorderCol;
    number[] mViscosityBorderRow;
    number[] mViscosityScratch;
    number[] mViscositySolution;
    int[] mViscosityPivots;
    number[] mThermalBorderCol;
    number[] mThermalBorderRow;
    number[] mThermalScratch;
    number[] mThermalSolution;
    int[] mThermalPivots;

    void allocateMatrix(ref number[][] matrix, int size)
    {
        matrix.length = size;
        foreach (row; 0 .. size) {
            matrix[row].length = size;
        }
    }
}

}

class ScaledTwoTemperatureTransProps : TwoTemperatureTransProps {
public:
    this(lua_State *L, string[] speciesNames)
    {
        super(L, speciesNames);
        mMuScale = getGlobalDoubleWithDefault(L, "mu_scale", 1.0);
        mKScale = getGlobalDoubleWithDefault(L, "k_scale", 1.0);
        mKModesScale = getGlobalDoubleWithDefault(L, "k_modes_scale", 1.0);
        mDiffusionScale = getGlobalDoubleWithDefault(L, "diffusion_scale", 1.0);
        requirePositive("mu_scale", mMuScale);
        requirePositive("k_scale", mKScale);
        requirePositive("k_modes_scale", mKModesScale);
        requirePositive("diffusion_scale", mDiffusionScale);
    }

    @nogc
    override void updateTransProps(ref GasState gs, GasModel gm)
    {
        super.updateTransProps(gs, gm);
        gs.mu *= mMuScale;
        gs.k *= mKScale;
        if (gs.k_modes.length > 0) gs.k_modes[0] *= mKModesScale;
    }

    @nogc
    override void binaryDiffusionCoefficients(ref const(GasState) gs, ref number[][] D)
    {
        super.binaryDiffusionCoefficients(gs, D);
        foreach (ref Di; D) {
            foreach (ref Dij; Di) {
                Dij *= mDiffusionScale;
            }
        }
    }

private:
    double mMuScale = 1.0;
    double mKScale = 1.0;
    double mKModesScale = 1.0;
    double mDiffusionScale = 1.0;

    double getGlobalDoubleWithDefault(lua_State *L, string key, double defaultValue)
    {
        double value = defaultValue;
        lua_getglobal(L, key.toStringz);
        if (lua_isnumber(L, -1)) value = lua_tonumber(L, -1);
        lua_pop(L, 1);
        return value;
    }

    void requirePositive(string key, double value)
    {
        if (value <= 0.0) {
            throw new Error(key ~ " must be positive for ScaledGuptaYos transport_model.");
        }
    }
}


unittest {
    import std.file;
    import gas.composite_gas;

    // Find the lua files in sample-data
    chdir("./sample-data");
    scope(exit) chdir("..");

    FloatingPointControl fpctrl;
    // Enable hardware exceptions for division by zero, overflow to infinity,
    // invalid operations, and uninitialized floating-point variables.
    // Copied from https://dlang.org/library/std/math/floating_point_control.html
    fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

    auto L = init_lua_State();
    GasModel gm = new CompositeGas("N2-2sp-2T.lua");
    doLuaFile(L, "N2-2sp-2T.lua");
    string[] speciesNames;
    getArrayOfStrings(L, "species", speciesNames);
    auto ttp = new TwoTemperatureTransProps(L, speciesNames);
    lua_close(L);
    auto gs = GasState(2, 1);
    gs.p = 1.0e5;
    gs.T = 2000.0;
    gs.T_modes[0] = 3000.0;
    gs.massf[0] = 0.8;
    gs.massf[1] = 0.2;

    ttp.updateTransProps(gs, gm);
    assert(approxEqualNumbers(to!number(6.5398632009870612e-05), gs.mu, 1.0e-2));
    assert(approxEqualNumbers(to!number(0.11258026890719305), gs.k, 1.0e-2));
    assert(approxEqualNumbers(to!number(0.018734503743319659), gs.k_modes[0], 1.0e-2));

    //import std.stdio;
    //writefln("mu= %.6e  k= %.6e  k_v= %.6e\n", gs.mu, gs.k, gs.k_modes[0]);

    gs.T_modes[0] = 300.0;
    ttp.updateTransProps(gs, gm);
    //writefln("mu = %.6e k = %.6e k_v = %.6e", gs.mu, gs.k, gs.k_modes[0]);
    assert(approxEqualNumbers(to!number(5.7677760411478017e-05), gs.k_modes[0], 1.0e-2));


    //auto L = init_lua_State();
    //GasModel gm = new CompositeGas("gm-air11-2T.lua");
    //doLuaFile(L, "gm-air11-2T.lua");
    //string[] speciesNames;
    //getArrayOfStrings(L, "species", speciesNames);
    //auto ttp = new TwoTemperatureTransProps(L, speciesNames);
    //lua_close(L);
    //auto gs = GasState(11, 1);
    //gs.p = 0.1*101.35e3;
    //gs.T = 1000.0;
    //gs.T_modes[0] = 3000.0;

    //gs.massf[0] =     7.455871e-01; // N2
    //gs.massf[1] =     2.543794e-01; // O2
    //gs.massf[2] =     0.000000e+00; // N
    //gs.massf[3] =     1.310344e-10; // O
    //gs.massf[4] =     3.355965e-05; // NO
    //gs.massf[5] =     0.000000e+00; // N2+
    //gs.massf[6] =     0.000000e+00; // O2+
    //gs.massf[7] =     0.000000e+00; // N+
    //gs.massf[8] =     0.000000e+00; // O+
    //gs.massf[9] =     0.000000e+00; // NO+
    //gs.massf[10] =     0.000000e+00; // e-


    //gm.update_thermo_from_pT(gs);
    //ttp.updateTransProps(gs, gm);
    //import std.stdio;
    //writefln("mu= %.6e  k= %.6e  k_v= %.6e\n", gs.mu, gs.k, gs.k_modes[0]);
}

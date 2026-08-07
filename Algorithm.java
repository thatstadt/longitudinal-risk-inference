import java.io.*;
import java.util.*;
import java.util.stream.*;

public class Algorithm {
    // -------------------------------------------------------------------------
    // Data generation
    // -------------------------------------------------------------------------
    private static final Random genRng = new Random(42);
    private static final double[] BMI_DRIFTS = {-0.2, -0.1, 0.0, 0.1, 0.2};
    private static final double[] setupFactors = new double[6];
    static {
        for (int i = 0; i < setupFactors.length; i++) setupFactors[i] = genRng.nextGaussian();
    }
    private static final int N               = 10_000_000;
    private static final int YEARS           = 25;
    private static final int BURNIN_YEARS    = 50;
    private static final int ITER            = 250;
    // Each bootstrap iteration resamples a FIXED subsample of the full
    // population (m-out-of-n bootstrap), rather than resampling all N
    // patients. With N now bumped to 10M, a full-N bootstrap would produce
    // extremely tight CIs (variance ~ 1/N); subsampling at a smaller, fixed
    // size estimates the sampling variability of a more modest study size
    // instead. The large N=10M population stabilizes the underlying "true"
    // simulated rate (law of large numbers); BOOT_SAMPLE_SIZE controls how
    // large a *real-world study* drawn from that true rate would be, so
    // this is the knob to shrink toward realistic clinical cohort sizes
    // (e.g. NHANES-scale, tens of thousands) for CIs comparable to GBD's.
    private static final int BOOT_SAMPLE_SIZE = 2_000_000;
    private static final int ID_LEN          = String.valueOf(N).length();

    private static final String[] GROUP_LABELS = {"0-14", "15-24", "25-44", "45-64", "65-80", "81+"};
    private static final int      NUM_GROUPS   = GROUP_LABELS.length;

    // US SSA 2019 period life table — blended male/female annual probability of death by age
    private static final double[] QX = {
        0.00535, 0.00036, 0.00023, 0.00018, 0.00014, // 0-4
        0.00013, 0.00012, 0.00011, 0.00010, 0.00009, // 5-9
        0.00009, 0.00010, 0.00012, 0.00016, 0.00022, // 10-14
        0.00027, 0.00033, 0.00039, 0.00045, 0.00051, // 15-19
        0.00057, 0.00062, 0.00067, 0.00071, 0.00073, // 20-24
        0.00074, 0.00076, 0.00079, 0.00083, 0.00088, // 25-29
        0.00093, 0.00100, 0.00107, 0.00116, 0.00126, // 30-34
        0.00137, 0.00149, 0.00163, 0.00178, 0.00196, // 35-39
        0.00215, 0.00237, 0.00262, 0.00289, 0.00318, // 40-44
        0.00350, 0.00386, 0.00426, 0.00471, 0.00520, // 45-49
        0.00575, 0.00635, 0.00701, 0.00771, 0.00844, // 50-54
        0.00923, 0.01006, 0.01098, 0.01200, 0.01313, // 55-59
        0.01438, 0.01574, 0.01720, 0.01875, 0.02041, // 60-64
        0.02221, 0.02417, 0.02631, 0.02864, 0.03117, // 65-69
        0.03391, 0.03684, 0.03997, 0.04331, 0.04688, // 70-74
        0.05073, 0.05490, 0.05942, 0.06437, 0.06977, // 75-79
        0.07567, 0.08212, 0.08914, 0.09674, 0.10492, // 80-84
        0.11367, 0.12297, 0.13278, 0.14306, 0.15378, // 85-89
        0.16489, 0.17634, 0.18809, 0.20009, 0.21227, // 90-94
        0.22459, 0.23699, 0.24942, 0.26182, 1.00000  // 95-99
    };
    private static final double RA_SMR = 1.1; // excess mortality for RA patients

    // Base names for output files — temp files are base + "_sim{i}.csv"
    private static final String[] FILE_BASES = {
        "prev_all", "inci_all", "prev_young", "inci_young",
        "prev_longitudinal_by_cohort", "inci_longitudinal_by_cohort"
    };


    private static int ageGroup(int age) {
        if (age < 15) return 0;
        if (age < 25) return 1;
        if (age < 45) return 2;
        if (age < 65) return 3;
        if (age < 81) return 4;
        return 5;
    }

    // -------------------------------------------------------------------------
    // Age-dependent baseline risk: f(age)
    //
    // Previously the model used a flat intercept (-9.34) plus a hand-tuned
    // logistic "a_beta" curve for age. That captured a rough rise-with-age
    // trend but could never reproduce the late-life *decline* in incidence
    // that GBD shows (incidence peaks ~65-69 then falls), because the
    // logistic just saturates at high age instead of turning over.
    //
    // Here we replace both the intercept and a_beta with a single smooth
    // function f(age), fit directly to GBD 2023 age-specific data:
    //
    //   1) For each 5-year GBD age band, take the target *annual incidence
    //      rate* I* and convert it to a target log-odds: logit(I*) =
    //      ln(I* / (1 - I*)). Incidence (not prevalence!) is the right
    //      calibration target here, because getRisk() returns an annual
    //      HAZARD -- the probability of onset in that one year -- and
    //      incidence is already expressed as a per-year rate. Prevalence,
    //      by contrast, is a *cumulative* result of decades of exposure to
    //      that hazard plus mortality/replacement, so calibrating f(age)
    //      directly against prevalence would (and on a first attempt,
    //      did) make the annual hazard itself far too large: e.g. forcing
    //      a ~1.25% chance of onset *every single year* at age 72 just
    //      because 1.25% of 72-year-olds carry the diagnosis explodes
    //      into ~20%+ cumulative prevalence after decades of exposure.
    //   2) Subtract the *expected* contribution of the non-age risk terms
    //      (genetics, sex, smoking) at that age. These population averages
    //      are known analytically from the model's own odds ratios and the
    //      generator's allele/sex/smoking rates (BMI is symmetric around
    //      the population mean and contributes ~0 in expectation). What's
    //      left is the age-only baseline log-odds required at that age.
    //   3) Fit a degree-5 least-squares polynomial through these
    //      (age midpoint, baseline log-odds) data points.
    //
    // The resulting polynomial naturally reproduces GBD's inverted-U shape:
    // low risk in youth, a peak in the mid-to-late 60s, and a decline beyond
    // that -- all driven by data rather than a hand-picked curve shape.
    // Prevalence then becomes an *emergent* property of the simulation
    // (cumulative hazard exposure + mortality + replacement), exactly as it
    // was under the old flat-intercept model -- compare the per-band
    // console/CSV prevalence output against GBD_PREV_PER_100K (kept below
    // for reference/validation only, not used in the fit) to sanity-check.
    //
    // CALIBRATION NOTE: step 2 ignores the small, run-specific constant
    // contributions of setupFactors[]/BMI noise, so simulated age-band
    // incidence will land close to but not exactly on the GBD targets. If
    // there's a roughly constant offset across bands, it can be absorbed
    // by nudging GBD_INCI_PER_100K uniformly and re-running -- the fit
    // itself does not need to change shape, only its overall level.
    // -------------------------------------------------------------------------

    // Age (years) at the center of each GBD 5-year band used for calibration.
    private static final double[] AGE_CALIB_MIDPOINTS = {
        7, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97
    };

    // GBD 2023 RA annual incidence rate, per 100,000 population, for the
    // band centered at the matching index in AGE_CALIB_MIDPOINTS. This is
    // the actual calibration target for f(age) -- see note above.
    private static final double[] GBD_INCI_PER_100K = {
        4.09,  // 0-14
        4.83,  // 15-19
        5.46,  // 20-24
        9.30,  // 25-29
        14.29, // 30-34
        20.40, // 35-39
        26.17, // 40-44
        31.44, // 45-49
        37.38, // 50-54
        44.16, // 55-59
        53.31, // 60-64
        62.64, // 65-69
        55.70, // 70-74
        34.31, // 75-79
        20.74, // 80-84
        14.79, // 85-89
        8.70,  // 90-94
        2.23   // 95+
    };

    // GBD 2023 RA prevalence, per 100,000 population -- NOT used in the fit.
    // Kept only so the simulation's emergent per-band prevalence output can
    // be compared against it as an independent validation check.
    private static final double[] GBD_PREV_PER_100K = {
        16.40,   // 0-14
        66.23,   // 15-19
        83.97,   // 20-24
        113.70,  // 25-29
        162.28,  // 30-34
        235.48,  // 35-39
        334.00,  // 40-44
        451.69,  // 45-49
        585.85,  // 50-54
        743.00,  // 55-59
        927.29,  // 60-64
        1125.12, // 65-69
        1245.33, // 70-74
        1220.64, // 75-79
        1095.49, // 80-84
        946.18,  // 85-89
        803.36,  // 90-94
        650.27   // 95+
    };

    // Coefficients of the calibrated degree-5 polynomial, lowest order first,
    // evaluated against (age / 100) for numerical stability. Computed once
    // at class load by fitting AGE_CALIB_MIDPOINTS -> calibration targets.
    private static final double[] AGE_POLY_COEFFS;

    static {
        int n = AGE_CALIB_MIDPOINTS.length;
        double[] scaledAge = new double[n];
        double[] target    = new double[n];

        // Allele-count levels with probabilities matching Patient() generator
        int[]    alleleCounts = {0, 1, 2};
        double[] alleleProbs  = {0.55, 0.37, 0.08};

        for (int i = 0; i < n; i++) {
            double age    = AGE_CALIB_MIDPOINTS[i];
            double iStar  = GBD_INCI_PER_100K[i] / 100_000.0;
            double targetLogit = Math.log(iStar / (1.0 - iStar));

            // We need to subtract E[exp(non-age log-odds)] -- in *log* space --
            // from the target logit to isolate the age baseline.
            //
            // The naive approach (subtract E[non-age betas]) is WRONG here:
            // for rare-event hazards, sigmoid(z) ≈ exp(z), so the population-
            // average hazard is exp(f_age) · E[exp(non-age terms)], NOT
            // exp(f_age + E[non-age terms]).  Those differ by Jensen's
            // inequality whenever per-patient terms vary across individuals --
            // the previous attempt that naively summed E[betas] produced targets
            // that were ~0.9 too low, causing a ~3.5x overshoot in the
            // resulting simulation.
            //
            // Instead we compute E[exp(non-age terms)] exactly by enumerating
            // all 12 (smoker × alleleCount × sex) combinations with their
            // joint probabilities -- analogous to a mixture-model integration.
            // For BMI (continuous Gaussian), E[exp(ln(b_eff)·Z)] for Z~N(0,1)
            // equals exp((ln(b_eff))^2 / 2) by the standard normal MGF.
            // The fixed setupFactors[] noise (seeded N(0,1)) is included
            // using the *actual* values already computed by the time this block
            // runs (they are declared and initialised above us in the file).
            double smokerRate  = age < 25 ? 0.05 : (age < 55 ? 0.13 : 0.10);
            double bEff        = 1.1 + 0.4 * Math.exp(-0.065 * Math.max(age - 15, 0));
            double bmiMgf      = Math.exp(Math.pow(Math.log(bEff), 2) / 2.0);

            double expExpNonAge = 0.0;
            for (int s = 0; s <= 1; s++) {
                double sProb  = s == 1 ? smokerRate : (1.0 - smokerRate);
                double sBeta  = s == 0 ? 0.0 : Math.log(2.5);
                for (int ai = 0; ai < alleleCounts.length; ai++) {
                    int    allele = alleleCounts[ai];
                    double aProb  = alleleProbs[ai];
                    double gOr    = allele == 2 ? 5.5 : (allele == 1 ? 2.0 : 1.0);
                    double gBeta  = Math.log(gOr);
                    double iOr    = 1.5 * s * (allele == 2 ? 1.0 : (allele == 1 ? 0.5 : 0.0));
                    double iBeta  = iOr == 0.0 ? 0.0 : Math.log(iOr);
                    for (int sx = 0; sx <= 1; sx++) {
                        double sxBeta = sx == 0 ? 0.0 : Math.log(2.5);
                        double betaSum = sBeta + gBeta + iBeta + sxBeta
                            + setupFactors[0] * 0.5  // smoking noise
                            + setupFactors[1] * 0.5  // genetics noise
                            + setupFactors[2] * 0.5  // interaction noise
                            + setupFactors[4] * 0.5; // sex noise
                        expExpNonAge += sProb * aProb * 0.5 * Math.exp(betaSum);
                    }
                }
            }

            // ln(E[exp(non-age)]) is the exact offset to remove from the
            // target logit so that f(age) isolates the age baseline correctly.
            double nonAgeE = Math.log(expExpNonAge * bmiMgf);

            scaledAge[i] = age / 100.0;
            target[i]    = targetLogit - nonAgeE;
        }

        AGE_POLY_COEFFS = fitPolynomial(scaledAge, target, 5);
    }

    /**
     * Fits a degree-{@code degree} polynomial y = c0 + c1*x + c2*x^2 + ...
     * to the given data points via ordinary least squares, solving the
     * normal equations (A^T A) c = A^T y with Gaussian elimination and
     * partial pivoting (A being the Vandermonde design matrix). Returns
     * the coefficients lowest-order first.
     */
    private static double[] fitPolynomial(double[] x, double[] y, int degree) {
        int n = x.length;
        int m = degree + 1;

        double[][] ata = new double[m][m];
        double[]   aty = new double[m];
        for (int i = 0; i < n; i++) {
            double[] powers = new double[2 * m - 1];
            powers[0] = 1.0;
            for (int p = 1; p < powers.length; p++) powers[p] = powers[p - 1] * x[i];
            for (int r = 0; r < m; r++) {
                aty[r] += powers[r] * y[i];
                for (int c = 0; c < m; c++) ata[r][c] += powers[r + c];
            }
        }

        for (int col = 0; col < m; col++) {
            int pivot = col;
            for (int r = col + 1; r < m; r++)
                if (Math.abs(ata[r][col]) > Math.abs(ata[pivot][col])) pivot = r;
            double[] tmpRow = ata[col]; ata[col] = ata[pivot]; ata[pivot] = tmpRow;
            double tmp = aty[col]; aty[col] = aty[pivot]; aty[pivot] = tmp;

            for (int r = col + 1; r < m; r++) {
                double factor = ata[r][col] / ata[col][col];
                for (int c = col; c < m; c++) ata[r][c] -= factor * ata[col][c];
                aty[r] -= factor * aty[col];
            }
        }

        double[] coeffs = new double[m];
        for (int r = m - 1; r >= 0; r--) {
            double sum = aty[r];
            for (int c = r + 1; c < m; c++) sum -= ata[r][c] * coeffs[c];
            coeffs[r] = sum / ata[r][r];
        }
        return coeffs;
    }

    /** Evaluates the calibrated age-baseline polynomial f(age) via Horner's method. */
    private static double agePoly(double age) {
        double xs = age / 100.0;
        double result = 0.0;
        for (int i = AGE_POLY_COEFFS.length - 1; i >= 0; i--)
            result = result * xs + AGE_POLY_COEFFS[i];
        return result;
    }

    private static class Patient {
        int age, smoker, alleleCount, sex, diagnosisYear, initialAgeGroup;
        double bmi, risk;
        boolean condition;

        Patient() {
            double chance = genRng.nextDouble();
            if      (chance < 0.18) age = genRng.nextInt(15);
            else if (chance < 0.94) age = genRng.nextInt(66) + 15;
            else                    age = genRng.nextInt(19) + 81;

            initialAgeGroup = ageGroup(age);
            diagnosisYear   = -1;
            condition       = false;

            chance = genRng.nextDouble();
            if      (age < 25) smoker = chance < 0.05 ? 1 : 0;
            else if (age < 55) smoker = chance < 0.13 ? 1 : 0;
            else               smoker = chance < 0.10 ? 1 : 0;

            chance      = genRng.nextDouble();
            alleleCount = chance < 0.08 ? 2 : (chance < 0.45 ? 1 : 0);
            sex         = genRng.nextDouble() < 0.5 ? 1 : 0;
            bmi         = Math.min(Math.max(genRng.nextGaussian() * 5.0 + 27.0, 15), 50);
            risk        = getRisk();
        }

        Patient(int age, double bmi, int smoker, int alleleCount, int sex) {
            this.age             = age;
            this.bmi             = bmi;
            this.smoker          = smoker;
            this.alleleCount     = alleleCount;
            this.sex             = sex;
            this.condition       = false;
            this.diagnosisYear   = -1;
            this.initialAgeGroup = ageGroup(age);
            this.risk            = getRisk();
        }

        double getRisk() {
            double s_beta   = smoker == 0 ? 0.0 : Math.log(2.5);
            double g_or     = alleleCount == 2 ? 5.5 : (alleleCount == 1 ? 2.0 : 1.0);
            double g_beta   = Math.log(g_or);
            double i_or     = 1.5 * smoker * (alleleCount == 2 ? 1.0 : (alleleCount == 1 ? 0.5 : 0.0));
            double i_beta   = i_or == 0.0 ? 0.0 : Math.log(i_or);
            double f_age    = agePoly(age); // calibrated age baseline, replaces intercept + a_beta
            double sx_beta  = sex == 0 ? 0.0 : Math.log(2.5);
            double bmiZ     = (bmi - 27.0) / 5.0;
            double b_eff    = 1.1 + 0.4 * Math.exp(-0.065 * Math.max(age - 15, 0));
            double b_beta   = Math.log(b_eff) * bmiZ;
            double z = f_age
                + s_beta  + setupFactors[0] * 0.5
                + g_beta  + setupFactors[1] * 0.5
                + i_beta  + setupFactors[2] * 0.5
                + sx_beta + setupFactors[4] * 0.5
                + b_beta;
            return 1.0 / (1.0 + Math.exp(-z));
        }

        void advanceOneYear(int year, Random rng, double bmiDrift) {
            age++;
            risk = getRisk();

            // Always consume the same two draws regardless of branch outcome
            // (common random numbers). Otherwise, the instant one patient's
            // onset outcome differs between two BMI-drift scenarios that
            // share a seed, every later rng call for every later patient in
            // that scenario desynchronizes from the others -- turning what
            // should be a clean, low-noise paired comparison across drift
            // scenarios into one dominated by spurious RNG-divergence noise.
            double onsetDraw = rng.nextDouble();
            double bmiDraw   = rng.nextGaussian();

            if (condition) return;
            if (onsetDraw < risk) {
                condition     = true;
                diagnosisYear = year;
                return;
            }
            bmi = Math.min(Math.max(bmi + bmiDraw * 0.3 + bmiDrift, 15), 50);
        }
    }

    private static void generateDataset(Patient[] patients) throws Exception {
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("dataset.csv")))) {
            pw.println("id,age,bmi,smoker,alleleCount,sex");
            for (int i = 0; i < patients.length; i++) {
                Patient p = patients[i];
                pw.printf("%0" + ID_LEN + "d,%d,%.6f,%d,%d,%d%n",
                    i + 1, p.age, p.bmi, p.smoker, p.alleleCount, p.sex);
            }
        }
        System.out.println("Dataset generated: dataset.csv");
    }

    private static Patient[] loadDataset() throws Exception {
        List<Patient> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader("dataset.csv"))) {
            br.readLine();
            String line;
            while ((line = br.readLine()) != null) {
                String[] p = line.split(",");
                list.add(new Patient(
                    Integer.parseInt(p[1]), Double.parseDouble(p[2]),
                    Integer.parseInt(p[3]), Integer.parseInt(p[4]), Integer.parseInt(p[5])
                ));
            }
        }
        System.out.println("Dataset loaded: dataset.csv (" + list.size() + " patients)");
        return list.toArray(new Patient[0]);
    }

    // -------------------------------------------------------------------------
    // Per-simulation runner (one thread per drift value)
    // -------------------------------------------------------------------------
    private static void runSimulation(
            double drift, int simIdx,
            int[] initAge, double[] initBmi, int[] initSmoker,
            int[] initAllele, int[] initSex) throws Exception {

        String driftHeader = "bmi_drift,year,avg,bmi_avg,delta,percent_increase";
        String longHeader  = "bmi_drift,year,group,avg,bmi_avg,nonbmi_avg,delta,percent_increase";

        String pA  = "prev_all_sim"                    + simIdx + ".csv";
        String iA  = "inci_all_sim"                    + simIdx + ".csv";
        String pY  = "prev_young_sim"                  + simIdx + ".csv";
        String iY  = "inci_young_sim"                  + simIdx + ".csv";
        String pL  = "prev_longitudinal_by_cohort_sim" + simIdx + ".csv";
        String iL  = "inci_longitudinal_by_cohort_sim" + simIdx + ".csv";

        try (PrintWriter prevAllPw  = new PrintWriter(new BufferedWriter(new FileWriter(pA)));
             PrintWriter inciAllPw  = new PrintWriter(new BufferedWriter(new FileWriter(iA)));
             PrintWriter prevYngPw  = new PrintWriter(new BufferedWriter(new FileWriter(pY)));
             PrintWriter inciYngPw  = new PrintWriter(new BufferedWriter(new FileWriter(iY)));
             PrintWriter prevLongPw = new PrintWriter(new BufferedWriter(new FileWriter(pL)));
             PrintWriter inciLongPw = new PrintWriter(new BufferedWriter(new FileWriter(iL)))) {

            prevAllPw.println(driftHeader);
            inciAllPw.println(driftHeader);
            prevYngPw.println(driftHeader);
            inciYngPw.println(driftHeader);
            prevLongPw.println(longHeader);
            inciLongPw.println(longHeader);

            // Restore cohort to initial state
            Patient[] patients = new Patient[N];
            for (int i = 0; i < N; i++)
                patients[i] = new Patient(initAge[i], initBmi[i], initSmoker[i], initAllele[i], initSex[i]);

            Random simRng     = new Random(99);
            Random replaceRng = new Random(77);

            // Burn-in: establish baseline prevalence, no drift, no recording
            for (int b = 0; b < BURNIN_YEARS; b++) {
                for (int i = 0; i < N; i++) patients[i].advanceOneYear(0, simRng, 0.0);
                for (int i = 0; i < N; i++) {
                    double qx = QX[Math.min(patients[i].age, 99)];
                    if (patients[i].condition) qx = Math.min(qx * RA_SMR, 1.0);
                    // Same common-random-numbers fix as advanceOneYear: always draw
                    // both qx-check and replacement-source index, regardless of
                    // whether a death occurs, to keep replaceRng synchronized
                    // across BMI-drift scenarios that share a seed.
                    boolean dies = simRng.nextDouble() < qx;
                    int src = replaceRng.nextInt(N);
                    if (dies) {
                        patients[i] = new Patient(initAge[src], initBmi[src], initSmoker[src], initAllele[src], initSex[src]);
                    }
                }
            }

            System.out.printf("%nSimulation: BMI drift = %+.1f%n", drift);
            System.out.printf("  Calibration targets -- inci(all): 0.02435%%  prev(all): 0.43911%%%n");

            for (int year = 1; year <= YEARS; year++) {

                final boolean[] preCondition = new boolean[N];
                for (int i = 0; i < N; i++) preCondition[i] = patients[i].condition;

                final double driftFinal = drift;
                for (int i = 0; i < N; i++) patients[i].advanceOneYear(year, simRng, driftFinal);

                final boolean[] cond      = new boolean[N];
                final boolean[] wasAtRisk = new boolean[N];
                final boolean[] newCase   = new boolean[N];
                final double[]  bmiA      = new double[N];
                final int[]     ageA      = new int[N];
                for (int i = 0; i < N; i++) {
                    cond[i]      = patients[i].condition;
                    wasAtRisk[i] = !preCondition[i];
                    newCase[i]   = patients[i].diagnosisYear == year;
                    bmiA[i]      = patients[i].bmi;
                    ageA[i]      = patients[i].age;
                }

                for (int i = 0; i < N; i++) {
                    double qx = QX[Math.min(patients[i].age, 99)];
                    if (patients[i].condition) qx = Math.min(qx * RA_SMR, 1.0);
                    // Same common-random-numbers fix as advanceOneYear: always draw
                    // both qx-check and replacement-source index, regardless of
                    // whether a death occurs, to keep replaceRng synchronized
                    // across BMI-drift scenarios that share a seed.
                    boolean dies = simRng.nextDouble() < qx;
                    int src = replaceRng.nextInt(N);
                    if (dies) {
                        // A population-wide BMI drift should affect everyone alive at
                        // a given simulated year, not just patients who have survived
                        // continuously since year 1 without ever dying and being
                        // replaced. Resetting replacements to the frozen pre-drift
                        // initBmi[] snapshot would dilute the drift signal exactly in
                        // the highest-mortality (oldest) age bands, where replacement
                        // turnover is highest -- so we shift the replacement's
                        // starting BMI by the drift accumulated since simulation
                        // start. At driftFinal == 0.0 (the reference scenario used
                        // for all GBD calibration comparisons) this offset is always
                        // exactly 0, so that scenario is completely unaffected.
                        double replacementBmi = Math.min(Math.max(initBmi[src] + driftFinal * year, 15), 50);
                        patients[i] = new Patient(initAge[src], replacementBmi, initSmoker[src], initAllele[src], initSex[src]);
                    }
                }

                double[][] prevAllRes  = new double[ITER][4];
                double[][] inciAllRes  = new double[ITER][4];
                double[][] prevYngRes  = new double[ITER][4];
                double[][] inciYngRes  = new double[ITER][4];
                double[][][] prevLongRes = new double[ITER][NUM_GROUPS][5];
                double[][][] inciLongRes = new double[ITER][NUM_GROUPS][5];

                final int yearFinal = year;

                IntStream.range(0, ITER).parallel().forEach(iter -> {
                    Random rng = new Random((long) yearFinal * 1000 + iter);

                    int pAllTot=0, pAllCond=0, pBmiTot=0, pBmiCond=0;
                    int pYngTot=0, pYngCond=0, pYngBmiTot=0, pYngBmiCond=0;
                    int iAllRisk=0, iAllNew=0, iBmiRisk=0, iBmiNew=0;
                    int iYngRisk=0, iYngNew=0, iYngBmiRisk=0, iYngBmiNew=0;
                    int[] gTot=new int[NUM_GROUPS],     gCond=new int[NUM_GROUPS];
                    int[] gBmiTot=new int[NUM_GROUPS],  gBmiCond=new int[NUM_GROUPS];
                    int[] gRisk=new int[NUM_GROUPS],    gNew=new int[NUM_GROUPS];
                    int[] gBmiRisk=new int[NUM_GROUPS], gBmiNew=new int[NUM_GROUPS];
                    int[] gNonBmiTot=new int[NUM_GROUPS],  gNonBmiCond=new int[NUM_GROUPS];
                    int[] gNonBmiRisk=new int[NUM_GROUPS], gNonBmiNew=new int[NUM_GROUPS];

                    for (int i = 0; i < BOOT_SAMPLE_SIZE; i++) {
                        int     idx    = rng.nextInt(N);
                        boolean c      = cond[idx];
                        boolean atRisk = wasAtRisk[idx];
                        boolean nc     = newCase[idx];
                        double  b      = bmiA[idx];
                        int     a      = ageA[idx];
                        boolean obese  = b >= 30.0;
                        boolean young  = a < 25;
                        int     g      = ageGroup(a);

                        pAllTot++; if (c) pAllCond++;
                        if (obese) { pBmiTot++; if (c) pBmiCond++; }

                        if (atRisk) { iAllRisk++; if (nc) iAllNew++; }
                        if (atRisk && obese) { iBmiRisk++; if (nc) iBmiNew++; }

                        if (young) {
                            pYngTot++; if (c) pYngCond++;
                            if (obese) { pYngBmiTot++; if (c) pYngBmiCond++; }
                        }
                        if (young && atRisk) { iYngRisk++; if (nc) iYngNew++; }
                        if (young && atRisk && obese) { iYngBmiRisk++; if (nc) iYngBmiNew++; }

                        gTot[g]++; if (c) gCond[g]++;
                        if (obese) { gBmiTot[g]++; if (c) gBmiCond[g]++; }
                        else       { gNonBmiTot[g]++; if (c) gNonBmiCond[g]++; }
                        if (atRisk) {
                            gRisk[g]++; if (nc) gNew[g]++;
                            if (obese) { gBmiRisk[g]++; if (nc) gBmiNew[g]++; }
                            else       { gNonBmiRisk[g]++; if (nc) gNonBmiNew[g]++; }
                        }
                    }

                    prevAllRes[iter] = stats((double) pAllCond / pAllTot,
                                              pBmiTot > 0 ? (double) pBmiCond / pBmiTot : 0);
                    inciAllRes[iter] = stats(iAllRisk > 0 ? (double) iAllNew / iAllRisk : 0,
                                              iBmiRisk > 0 ? (double) iBmiNew / iBmiRisk : 0);
                    prevYngRes[iter] = stats(pYngTot    > 0 ? (double) pYngCond    / pYngTot    : 0,
                                              pYngBmiTot > 0 ? (double) pYngBmiCond / pYngBmiTot : 0);
                    inciYngRes[iter] = stats(iYngRisk    > 0 ? (double) iYngNew    / iYngRisk    : 0,
                                              iYngBmiRisk > 0 ? (double) iYngBmiNew / iYngBmiRisk : 0);
                    for (int g = 0; g < NUM_GROUPS; g++) {
                        prevLongRes[iter][g] = stats3(gTot[g]    > 0 ? (double) gCond[g]    / gTot[g]    : 0,
                                                       gBmiTot[g] > 0 ? (double) gBmiCond[g] / gBmiTot[g] : 0,
                                                       gNonBmiTot[g] > 0 ? (double) gNonBmiCond[g] / gNonBmiTot[g] : 0);
                        inciLongRes[iter][g] = stats3(gRisk[g]    > 0 ? (double) gNew[g]    / gRisk[g]    : 0,
                                                       gBmiRisk[g] > 0 ? (double) gBmiNew[g] / gBmiRisk[g] : 0,
                                                       gNonBmiRisk[g] > 0 ? (double) gNonBmiNew[g] / gNonBmiRisk[g] : 0);
                    }
                });

                StringBuilder sb = new StringBuilder();
                for (int it = 0; it < ITER; it++) {
                    writeRow(sb, prevAllPw,  drift, year, null,             prevAllRes[it]);
                    writeRow(sb, inciAllPw,  drift, year, null,             inciAllRes[it]);
                    writeRow(sb, prevYngPw,  drift, year, null,             prevYngRes[it]);
                    writeRow(sb, inciYngPw,  drift, year, null,             inciYngRes[it]);
                    for (int g = 0; g < NUM_GROUPS; g++) {
                        writeRow(sb, prevLongPw, drift, year, GROUP_LABELS[g], prevLongRes[it][g]);
                        writeRow(sb, inciLongPw, drift, year, GROUP_LABELS[g], inciLongRes[it][g]);
                    }
                }

                double[] mInciAll  = new double[ITER];
                double[] mInciYng  = new double[ITER];
                double[] mPrevAll  = new double[ITER];
                double[] mPrevYng  = new double[ITER];
                double[][] mInciLong = new double[NUM_GROUPS][ITER];
                double[][] mPrevLong = new double[NUM_GROUPS][ITER];
                for (int it = 0; it < ITER; it++) {
                    mInciAll[it] = inciAllRes[it][0];
                    mInciYng[it] = inciYngRes[it][0];
                    mPrevAll[it] = prevAllRes[it][0];
                    mPrevYng[it] = prevYngRes[it][0];
                    for (int g = 0; g < NUM_GROUPS; g++) {
                        mInciLong[g][it] = inciLongRes[it][g][0];
                        mPrevLong[g][it] = prevLongRes[it][g][0];
                    }
                }
                Arrays.sort(mInciAll); Arrays.sort(mInciYng);
                Arrays.sort(mPrevAll); Arrays.sort(mPrevYng);
                for (int g = 0; g < NUM_GROUPS; g++) {
                    Arrays.sort(mInciLong[g]);
                    Arrays.sort(mPrevLong[g]);
                }

                int lo = (int)(0.025 * ITER), mid = (int)(0.500 * ITER), hi = (int)(0.975 * ITER);
                int yng = 1;
                System.out.printf("  Year %2d  [inci(all): %.5f%% vs 0.02435%%  |  prev(all): %.4f%% vs 0.43911%%]%n",
                    year, mInciAll[mid], mPrevAll[mid]);
                System.out.printf("    inci: all=%.4f[%.4f,%.4f]  young=%.4f[%.4f,%.4f]  long(15-24)=%.4f[%.4f,%.4f]%n",
                    mInciAll[mid], mInciAll[lo], mInciAll[hi],
                    mInciYng[mid], mInciYng[lo], mInciYng[hi],
                    mInciLong[yng][mid], mInciLong[yng][lo], mInciLong[yng][hi]);
                System.out.printf("    prev: all=%.4f[%.4f,%.4f]  young=%.4f[%.4f,%.4f]  long(15-24)=%.4f[%.4f,%.4f]%n",
                    mPrevAll[mid], mPrevAll[lo], mPrevAll[hi],
                    mPrevYng[mid], mPrevYng[lo], mPrevYng[hi],
                    mPrevLong[yng][mid], mPrevLong[yng][lo], mPrevLong[yng][hi]);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Merge temp files into final output files, then delete temps
    // -------------------------------------------------------------------------
    private static void mergeAndCleanup() throws Exception {
        String driftHeader = "bmi_drift,year,avg,bmi_avg,delta,percent_increase";
        String longHeader  = "bmi_drift,year,group,avg,bmi_avg,nonbmi_avg,delta,percent_increase";

        for (String base : FILE_BASES) {
            String header = base.contains("longitudinal") ? longHeader : driftHeader;
            try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(base + ".csv")))) {
                out.println(header);
                for (int si = 0; si < BMI_DRIFTS.length; si++) {
                    File temp = new File(base + "_sim" + si + ".csv");
                    try (BufferedReader br = new BufferedReader(new FileReader(temp))) {
                        br.readLine(); // skip temp header
                        String line;
                        while ((line = br.readLine()) != null) out.println(line);
                    }
                    temp.delete();
                }
            }
            System.out.println("Merged: " + base + ".csv");
        }
    }

    // -------------------------------------------------------------------------
    // Main
    // -------------------------------------------------------------------------
    public static void main(String[] args) throws Exception {
        // Generate patients once — all simulations share the same initial cohort
        Patient[] seed = new Patient[N];
        for (int i = 0; i < N; i++) seed[i] = new Patient();
        generateDataset(seed); // Comment out to skip generation and load existing dataset.csv instead
        // seed = loadDataset(); // Uncomment to load from existing dataset.csv

        // Snapshot initial state so each simulation can restore it
        final int[]    initAge    = new int[N];
        final double[] initBmi    = new double[N];
        final int[]    initSmoker = new int[N];
        final int[]    initAllele = new int[N];
        final int[]    initSex    = new int[N];
        for (int i = 0; i < N; i++) {
            initAge[i]    = seed[i].age;
            initBmi[i]    = seed[i].bmi;
            initSmoker[i] = seed[i].smoker;
            initAllele[i] = seed[i].alleleCount;
            initSex[i]    = seed[i].sex;
        }
        seed = null; // allow GC — no longer needed

        // Run drift simulations sequentially — inner bootstrap already uses all cores
        for (int si = 0; si < BMI_DRIFTS.length; si++)
            runSimulation(BMI_DRIFTS[si], si, initAge, initBmi, initSmoker, initAllele, initSex);

        // Merge temp files into final CSVs and delete temps
        System.out.println("\nMerging output files...");
        mergeAndCleanup();

        System.out.println("Program finished.");
    }

    private static double[] stats(double avg, double bmiAvg) {
        double delta = bmiAvg - avg;
        double pct   = avg > 0 ? (bmiAvg / avg - 1.0) * 100 : 0;
        return new double[]{ avg * 100, bmiAvg * 100, delta * 100, pct };
    }

    // Like stats(), but also reports the non-obese-only rate, and computes
    // delta/percent_increase as obese vs. non-obese (a true two-group
    // comparison) rather than obese vs. the whole pooled group (which
    // dilutes the contrast, since the whole group already contains the
    // obese patients).
    private static double[] stats3(double avg, double bmiAvg, double nonBmiAvg) {
        double delta = bmiAvg - nonBmiAvg;
        double pct   = nonBmiAvg > 0 ? (bmiAvg / nonBmiAvg - 1.0) * 100 : 0;
        return new double[]{ avg * 100, bmiAvg * 100, nonBmiAvg * 100, delta * 100, pct };
    }

    private static void writeRow(StringBuilder sb, PrintWriter pw, double drift, int year, String group, double[] row) {
        sb.setLength(0);
        sb.append(drift);
        sb.append(',');
        sb.append(year);
        if (group != null) { sb.append(','); sb.append(group); }
        for (double v : row) { sb.append(','); sb.append(v); }
        pw.println(sb);
    }
}

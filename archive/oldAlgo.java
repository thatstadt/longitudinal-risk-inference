//Algorithm.java
import java.io.*;
import java.util.*;

public class oldAlgo {
    private static final Random rng = new Random(42);
    private static final HashMap<String, Double> setupFactors = new HashMap<>();
    static {
        setupFactors.put("smoking", rng.nextGaussian());
        setupFactors.put("genetics", rng.nextGaussian());
        setupFactors.put("interaction", rng.nextGaussian());
        setupFactors.put("age", rng.nextGaussian());
        setupFactors.put("sex", rng.nextGaussian());
        setupFactors.put("bmi", rng.nextGaussian());
    }
    private static int n = 100_000;
    private static int n_25 = 0;

    public static void main(String[] args) throws Exception {
        
        //Create data
        ArrayList<Patient> data = new ArrayList<>();
        for (int i = 1; i <= n; i++) {
            Patient p = new Patient(i);
            if (p.age < 25) n_25++;
            data.add(p);
        }
        exportPatients(data, "dataset.csv");
        data.clear();

        //Import data
        FileReader fr = new FileReader(new File("dataset.csv"));
        Scanner sc = new Scanner(fr);
        sc.nextLine(); //Skip header
        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            String[] parts = line.split(",");
            Patient p = new Patient(); //Dummy patient
            p.id = parts[0];
            p.age = Integer.parseInt(parts[1]);
            p.bmi = Double.parseDouble(parts[2]);
            p.risk = Double.parseDouble(parts[3]);
            p.condition = Boolean.parseBoolean(parts[4]);
            data.add(p);
        }
        sc.close();

        double[] lineData = new double[4];
        for (int ui = 0; ui < 100; ui++) {
            ArrayList<Patient> subset = new ArrayList<>();
            for (int i = 1; i <= n; i++) {
                int chance = rng.nextInt(n);
                subset.add(data.get(chance));
            }

            //Average, total
            double avgRisk = 0.0;
            double delta = 0.0;
            double pIncrease = 1.0;
            int total = 0;

            for (Patient p : subset) {
                avgRisk += p.condition ? 1 : 0;
                total++;
            }

            avgRisk /= total;
            delta -= avgRisk;
            pIncrease /= avgRisk;
            lineData[0] = avgRisk * 100;
            avgRisk = 0.0;
            total = 0;

            //Average, obesity-specific
            for (Patient p : subset) {
                if (p.bmi >= 30.0) {
                    avgRisk += p.condition ? 1 : 0;
                    total++;
                }
            }
            avgRisk /= total;
            delta += avgRisk;
            pIncrease *= avgRisk;
            pIncrease -= 1.0;
            lineData[1] = avgRisk * 100;
            avgRisk = 0.0;
            total = 0;
            lineData[2] = delta * 100;
            delta = 0.0;
            lineData[3] = pIncrease * 100;
            pIncrease = 1.0;

            //Export to file
            exportData(lineData, "all_ages.csv");

            //Average, young-onset specific
            for (Patient p : subset) {
                if (p.age < 25) {
                    avgRisk += p.condition ? 1 : 0;
                    total++;
                }
            }
            avgRisk /= total;
            delta -= avgRisk;
            pIncrease /= avgRisk;
            lineData[0] = avgRisk * 100;
            avgRisk = 0.0;
            total = 0;

            //Average, young-onset obesity specific
            for (Patient p : subset) {
                if (p.bmi >= 30 && p.age < 25) {
                    avgRisk += p.condition ? 1 : 0;
                    total++;
                }
            }
            avgRisk /= total;
            delta += avgRisk;
            pIncrease *= avgRisk;
            pIncrease -= 1.0;
            lineData[1] = avgRisk * 100;
            avgRisk = 0.0;
            total = 0;
            lineData[2] = delta * 100;
            delta = 0.0;
            lineData[3] = pIncrease * 100;
            pIncrease = 1.0;

            //Export to file
            exportData(lineData, "ag_1.csv");
        }

        fr = new FileReader(new File("ag_1.csv"));
        sc = new Scanner(fr);
        sc.nextLine();

        ArrayList<Double> values = new ArrayList<>();

        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            String[] parts = line.split(",");
            values.add(Double.parseDouble(parts[0]));
        }
        sc.close();

        Collections.sort(values);

        int size = values.size();

        double lower = values.get((int)(0.025 * size));
        double median = values.get((int)(0.5 * size));
        double upper = values.get((int)(0.975 * size));

        System.out.println("CI: median = " + median + ", lower = " + lower + ", upper = " + upper);
        System.out.print("Program finished.");
    }
    
    private static class Patient {
        //Patient attributes
        private String id;
        private int smoker; //1 if smoker, 0 if not
        private int age;
        private int alleleCount; //0, 1, or 2
        private int sex; //0 if male, 1 if female
        private double bmi;
        private double risk; //Probabilistic risk factor
        private boolean condition; //true if RA-positive, false if not
        private HashMap<String, Double[]> vars;
        private int diagnosisYear;

        public Patient() {}

        public Patient(int i) {
            //Constructs patient ID
            StringBuilder sb = new StringBuilder();
            while (sb.length() < String.valueOf(n).length() - String.valueOf(i).length()) {
                sb.append(0);
            }
            sb.append(i);
            this.id = sb.toString();

            //Age function
            double chance = rng.nextDouble();
            if (chance < 0.18) {
                this.age = rng.nextInt(15); //Age 0-14
            } else if (chance < 0.94) {
                this.age = rng.nextInt(66) + 15; //Age 15-80
            } else {
                this.age = rng.nextInt(19) + 81; //Age 81-99
            }

            //Diagnosis year setup
            this.diagnosisYear = -1;

            //Smoker atribute
            chance = rng.nextDouble();
            if (this.age < 25) this.smoker = chance < 0.05 ? 1 : 0;
            else if (this.age < 55) this.smoker = chance < 0.13 ? 1 : 0;
            else this.smoker = chance < 0.1 ? 1 : 0;

            //Allele count
            chance = rng.nextDouble();
            this.alleleCount = chance < 0.08 ? 2 : (chance < 0.45 ? 1 : 0);

            //Sex determination
            this.sex = rng.nextDouble() < 0.5 ? 1 : 0;

            //BMI setup
            this.bmi = Math.min(Math.max(rng.nextGaussian() * 5.0 + 27.0, 15), 50); //Mean 27, SD 5

            //Initial risk and condition status
            this.condition = false;
            for (int year = 0; year < 5; year++) {
                advanceOneYear(year);
            }
        }

        //Initializes the HashMap for variables
        public HashMap<String, Double[]> initializeVars() {
            HashMap<String, Double[]> init = new HashMap<>();

            //Smoking
            init.put("smoking", smoking());
            //init.put("smoking", new Double[]{0.0, 0.0});

            //Genetics
            init.put("genetics", genetics());
            //init.put("genetics", new Double[]{0.0, 0.0});

            //Smoking-genetics interaction
            init.put("interaction", interaction());
            //init.put("interaction", new Double[]{0.0, 0.0});

            //Age
            init.put("age", age());
            //init.put("age", new Double[]{0.0, 0.0});

            //Sex
            init.put("sex", sex());
            //init.put("sex", new Double[]{0.0, 0.0});

            //BMI
            init.put("bmi", bmi());
            //init.put("bmi", new Double[]{0.0, 0.0});

            return init;
        }

        private static double sigmoid(double z) {
            return 1.0 / (1.0 + Math.exp(-z));
        }

        //Smoking OR, std. deviation
        private Double[] smoking() {
            double s_or = 2.5 * this.smoker;
            double s_beta = Math.log(s_or == 0 ? 1 : s_or);
            double s_std = 0.5;
            Double[] s = {s_beta, s_std};
            return s;
        }

        //Genetics OR, std. deviation
        private Double[] genetics() {
            double g_or = this.alleleCount == 2 ? 5.5 : (this.alleleCount == 1 ? 2.0 : 1);
            double g_beta = Math.log(g_or == 0 ? 1 : g_or);
            double g_std = 0.5;
            Double[] g = {g_beta, g_std};
            return g;
        }

        //Smoking-genetics interaction OR, std. deviation
        private Double[] interaction() {
            double i_or = 1.5 * this.smoker * (this.alleleCount == 2 ? 1 : (this.alleleCount == 1 ? 0.5 : 0));
            double i_beta = Math.log(i_or == 0 ? 1 : i_or);
            double i_std = 0.5;
            Double[] i = {i_beta, i_std};
            return i;
        }

        //Age OR, std. deviation - WIP
        private Double[] age() {
            double minAge = 15.0;
            double maxAge = 80.0;
            double ageRange = maxAge - minAge;

            //Normalize age from 0 to 1
            double fraction = (this.age - minAge) / ageRange;
            fraction = Math.max(0.0, Math.min(fraction, 1.0));

            //Logistic curve parameters for young-age protection
            double steepness = 6.0; //Controls how quickly risk increases with age
            double midpoint = 0.5; //Midpoint of the curve (age ~45)

            double minBeta = -3.0; //Strong negative log-odds at young age <-----
            double maxBeta = 0.1; //Mild positive log-odds at old age <-----

            //Sigmoid for smooth transition
            double logistic = sigmoid(steepness * (fraction - midpoint));

            double a_beta = minBeta + logistic * (maxBeta - minBeta);
            double a_std = 0.0;
            Double[] a = {a_beta, a_std};
            return a;
        }

        //Sex OR, std. deviation
        private Double[] sex() {
            double s_or = 2.5 * this.sex;
            double s_beta = Math.log(s_or == 0 ? 1 : s_or);
            double s_std = 0.5;
            Double[] s = {s_beta, s_std};
            return s;
        }

        //BMI OR, std. deviation
        private Double[] bmi() {
            double k = -0.065;
            double bmiMean = 27.0;
            double bmiStd = 5.0;
            double bmiZ = (this.bmi - bmiMean) / bmiStd;

            double b_max = 1.26;
            double b_eff = Math.max(1.0 + (b_max - 1.0) * Math.exp(k * Math.max(this.age - 15, 0)), 1.01);
            
            double b_beta = Math.log(b_eff) * bmiZ;
            double b_std = 0.0;
            Double[] b = {b_beta, b_std};
            return b;
        }

        //Updates BMI each year with some random fluctuation
        private void updateBMI() {
            double yearlyChange = rng.nextGaussian() * 0.3 + 0.1;
            bmi = Math.min(Math.max(bmi + yearlyChange, 15), 50);
        }

        //Returns the betas of vars
        public double getRisk() {
            this.vars = initializeVars(); //Initialize variables

            double z = -6.0; //Intercept - find a real-world based value (Test training splits)
            for (String key : vars.keySet()) {
                Double[] v = vars.get(key);
                z += v[0] + setupFactors.get(key) * v[1];
            }
            return sigmoid(z);
        }

        //Advances patient by one year, updating age, risk, and condition status
        public void advanceOneYear(int year) {
            this.age++;
            this.risk = getRisk();
            
            if (this.condition) return;
            if (rng.nextDouble() < this.risk) {
                this.condition = true;
                this.diagnosisYear = year;
                return;
            }

            updateBMI();
        }
    }

    public static void exportPatients(ArrayList<Patient> patients, String filename) throws Exception {
        try (PrintWriter pw = new PrintWriter(filename)) {
            pw.println("id,age,bmi,risk,condition,diagnosisYear");
            for (Patient p : patients) {
                pw.println(p.id + "," + p.age + "," + p.bmi + "," + p.risk + "," + p.condition + "," + p.diagnosisYear);
            }
            pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void exportData(double[] lineData, String filename) throws Exception {
        File file = new File(filename);
        boolean header = file.exists();

        try (PrintWriter pw = new PrintWriter(new FileWriter(filename, true))) {
            if (!header || file.length() == 0) pw.println("avg,bmi_avg,delta,percent_increase"); //Header
            for (int i = 0; i < lineData.length; i++) {
                pw.print(lineData[i]);
                if (i != lineData.length - 1) pw.print(",");
            }
            pw.println();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

package de.unijena.cheminf.fragment.fingerprint.performanceTest;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
          //  PerformanceTest tmpApplication = new PerformanceTest("Fragments_Ertl_algorithm_200k_COCONUT.csv","Items_Ertl_algorithm_200k_COCONUT.csv","200000");
       // PerformanceTest tmpApplication = new PerformanceTest("Fragments_Ertl_algorithm_200k_COCONUT.csv","Items_Ertl_algorithm_200k_COCONUT.csv","200000");
            PerformanceTest tmpApplication = new PerformanceTest(args[0], args[1], args[2]);

    }
}

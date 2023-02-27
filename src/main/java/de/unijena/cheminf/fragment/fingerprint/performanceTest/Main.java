/*
 * MIT License
 *
 * Copyright (c) 2023 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.fragment.fingerprint.performanceTest;

/**
 * Main class for starting the PerformanceTest application.
 *
 * @author Betuel Sevindik
 */
public class Main {
    /**
     * Starts the application.Command line arguments must be the name of CSV-files to read.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args)  {
        try {
           // PerformanceTest tmpApplication = new PerformanceTest("Fragments_Ertl_algorithm_200k_COCONUT.csv", "Items_Ertl_algorithm_200k_COCONUT.csv","20000");
            PerformanceTest tmpApplication = new PerformanceTest(args[0], args[1], args[2]);
        } catch (Exception anException) {
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }
}

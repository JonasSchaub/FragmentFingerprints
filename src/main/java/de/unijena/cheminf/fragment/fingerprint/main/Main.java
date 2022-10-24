package de.unijena.cheminf.fragment.fingerprint.main;

import de.unijena.cheminf.fragment.fingerprint.BitFragmentFingerprint;
import de.unijena.cheminf.fragment.fingerprint.CountFragmentFingerprint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

public class Main {

    public static void main(String[] args) throws IOException {


        BitFragmentFingerprint bitFragmentFingerprint = new BitFragmentFingerprint();
        CountFragmentFingerprint countFragmentFingerprint = new CountFragmentFingerprint();

        BufferedReader br2 = new BufferedReader(new FileReader("C:\\Users\\betue\\Desktop\\Items_Ertl_algorithm_test.csv"));
        String line2 = "";
        String DELIMITER2 = ";";
        List<String> records2 = new ArrayList<>();
        List<List<String>> records3 = new ArrayList<>();
        HashMap<String, String> maps = new HashMap<>();
        String line3;
        while ((line3 = br2.readLine()) != null) {
            String[] values = line3.split(DELIMITER2);
            records3.add(Arrays.asList(values));

          //  System.out.println(records3.get(0) + "indexNull");
        }
        List<String> str = null;
        for(int z = 1; z< records3.size(); z++) {
            str = records3.get(z);
            //System.out.println(str + "---einezlene Listen");
            HashMap<String, String> allMaps = new HashMap<>();
            for(int i = 0; i<str.size(); i++) {
                if (i % 2 == 0) {
                    allMaps.put(str.get(i), str.get(i + 1));
                }
            }
           // System.out.println(allMaps +"------map");
            ArrayList g = bitFragmentFingerprint.generateFragmentFingerprint((ArrayList<String>)getFragmentList(new File("C:\\Users\\betue\\Desktop\\Fragments_Ertl_algorithm.csv")), allMaps);
          //  System.out.println( g.get(g.keySet().toArray()[0]));
            System.out.println(bitFragmentFingerprint.getNumberPositiveBits(g) + "neue Methode");
            int[] arr = (int[]) g.get(0);
            System.out.println(java.util.Arrays.toString(arr) +"--"+str.get(0));
            //System.out.println(java.util.Arrays.toString(g) +"--"+str.get(0));
           // int[] q = countFragmentFingerprint.generateFragmentFingerprint((ArrayList<String>) records, allMaps);
           // System.out.println(java.util.Arrays.toString(q)+ "---"+str.get(0));
        }

    }

    private static ArrayList<String> getFragmentList(File aFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(aFile));
        String line = "";
        String DELIMITER = ",";
        ArrayList<String> records = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            String[] values = line.split(DELIMITER);
            records.add(values[0]);
        }
        records.remove(0);
        return  records;
        // System.out.println(records+ "----fragment list");
    }
}

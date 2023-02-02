package TRM;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;

/**
 * @author lmh
 * @since 2022/11/2
 **/
public class Main {
    public static void main(String[] args) throws IOException {
//        String dataset = "10k100E32L";
//        String dataset = "20k100E32L";
//        String dataset = "ASL_BU_1";
        String dataset = "ASL_BU_2";
//        String dataset = "CONTEXT";
//        String dataset = "HEPATITIS";
        double[] us = new double[]{0.0017, 0.00175, 0.0018, 0.00185, 0.0019, 0.00195, 0.002};
        double utilityratio = 35;
        double confidence = 0.6;
        String datasetFile = "intervalData/" + dataset + ".txt";
        String input = fileToPath(datasetFile);
        String output = "rule in " + dataset + " thd " + utilityratio + ".txt";

        // This parameter let the user specify how many sequences from the input file should be used.
        // For example, it could be used to read only the first 1000 sequences of an input file
        int maximumSequenceCount = Integer.MAX_VALUE;

        //  THESE ARE ADDITIONAL PARAMETERS
        //   THE FIRST PARAMETER IS A CONSTRAINT ON THE MAXIMUM NUMBER OF INTERVALS IN THE LEFT SIDE OF RULES
        // For example, we don't want to find rules with more than 9 intervals in their left side
        int maxAntecedentSize = 9;  // 9

        //   THE SECOND PARAMETER IS A CONSTRAINT ON THE MAXIMUM NUMBER OF INTERVALS IN THE RIGHT SIDE OF RULES
        // For example, we don't want to find rules with more than 9 intervals in their right side
        int maxConsequentSize = 9;  // 9

        Algo algo = new Algo();
//        Algo_nmc algo = new Algo_nmc();
        algo.runAlgo(input, output, utilityratio, confidence, maxAntecedentSize, maxConsequentSize, maximumSequenceCount);
    }

    public static String fileToPath(String filename) throws UnsupportedEncodingException {
        URL url = Algo.class.getResource(filename);
        return java.net.URLDecoder.decode(((URL) url).getPath(),"UTF-8");
    }
}

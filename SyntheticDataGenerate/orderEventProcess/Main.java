package orderEventProcess;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URL;

/**
 * @Copyright (C), HITSZ
 * @Author Maohua Lv
 * @Date Created in 15:34 2022/11/8.
 * @Version 1.0
 * @Description
 */
public class Main {
    public static void main(String[] args) throws IOException {
        int sequenceCount = 40000;

        int maxDistinctItems = 100;

        int meanIntervalCountBySequence = 32;

        String output = sequenceCount / 1000 + "k" + maxDistinctItems + "E" + meanIntervalCountBySequence + "L" + ".txt";
        IntervalSequenceGenerator algo = new IntervalSequenceGenerator();
        algo.generateDatabase(sequenceCount, maxDistinctItems, meanIntervalCountBySequence, output);
    }

    public static String fileToPath(String filename) throws UnsupportedEncodingException {
        URL url = Algo.class.getResource(filename);
        return java.net.URLDecoder.decode(((URL) url).getPath(),"UTF-8");
    }
}

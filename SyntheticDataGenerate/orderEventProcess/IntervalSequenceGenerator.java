package orderEventProcess;

/* This file is copyright (c) 2008-2012 Philippe Fournier-Viger
 *
 * This file is part of the SPMF DATA MINING SOFTWARE
 * (http://www.philippe-fournier-viger.com/spmf).
 *
 * SPMF is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along with
 * SPMF. If not, see <http://www.gnu.org/licenses/>.
 */

import orderEventProcess.Interval;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * This class is a random sequence database generator such that
 * the user provides some parameters and this class generate a sequence database
 * that is written to the disk.
 *
 * @author Philippe Fournier-Viger
 */
public class IntervalSequenceGenerator {

	// a random number generator
	private static Random random = new Random(20003);

	/**
	 * This method randomly generates a sequence database according to parameters provided.
	 * @param sequenceCount the number of sequences required
	 * @param maxDistinctItems the maximum number of distinct items
	 * @param meanIntervalCountBySequence the average number of intervals by sequence
	 * @param intervalDataFile the file path for writting the generated database
	 * @throws IOException
	 */
	public void generateDatabase(int sequenceCount, int maxDistinctItems, int meanIntervalCountBySequence, String intervalDataFile) throws IOException {

		// We create a BufferedWriter to write the database to disk
		BufferedWriter writer = new BufferedWriter(new FileWriter(intervalDataFile));

		int k = 1;
		// For the number of sequences to be generated
		for (int i = 1; i <= sequenceCount; i++) {
			HashMap<Integer, Integer> evtTime = new HashMap<>();
			if (i % 1000 == 0) System.out.println("complete generate " + k++ + "k sequences");

			//start and finish time
			int starttime = 0;
			int endtime = 0;
			int gap = 13;
//			while (gap < 0 || gap > 30) {
//				gap = random.nextInt();
//			}
			//last interval event
//			int lastEvent = 0;
			//last event's start time
			int laststarttime = 1;

			//interval's number。variance 9，mean = meanIntervalCountBySequence
			double intervalNum1 = Math.abs(Math.sqrt(9) * random.nextGaussian() + meanIntervalCountBySequence);
			int intervalNum = new Double(intervalNum1).intValue();

			// for the number of intervals to be generated
			for (int j = 0; j < intervalNum; j++) {
				/**start time**/
				starttime = random.nextInt();
				while(starttime < laststarttime || starttime > endtime + gap) {
					starttime = random.nextInt(10000) + 1;
					System.out.println("seq: " + i + " old start time: " + laststarttime + " new start time: " + starttime);
				}
				if(laststarttime > starttime) starttime = laststarttime + 1; //to avoid get a smaller start time of the generated event
				laststarttime = starttime;

				/**finish time**/
				//duration. nomal distribution：variance，mean
				double mean = (random.nextInt(3) + 1) * 100;
				double durationDouble = Math.abs(Math.sqrt(1000) * random.nextGaussian() + mean);
				int duration = new Double(durationDouble).intValue();
				if(duration == 0) duration++;
				endtime = starttime + duration;

				/**utility**/
				//utility，log nomal	duration
				double utilityNormal = Math.sqrt(8) * random.nextGaussian() + Math.exp(Math.abs((1.0 / (Math.sqrt(2 *
						Math.PI * 8))) * Math.exp(-((duration - 8) * (duration - 8)) / (16.))));
				double utilityLogNormal = Math.exp(utilityNormal);
//				if(utilityLogNormal >= 80d) utilityLogNormal = Math.exp(Math.log(utilityNormal) + Math.E);
				int utility = (int)Math.ceil(utilityLogNormal);

				/**event**/
				int event = random.nextInt(maxDistinctItems) + 1;
				while(evtTime.get(event) != null && evtTime.get(event) >= endtime){ 
					event = random.nextInt(maxDistinctItems);
				}

				evtTime.put(event, Math.max(evtTime.getOrDefault(event, 0), endtime));

				StringBuilder buffer = new StringBuilder();
				buffer.append(i + " " + event + " " + starttime + " " + endtime + " " + utility);
				//buffer.append(i + "," + event + "," + starttime + "," + endtime);
				writer.write(buffer.toString());
				writer.newLine();
			}
		}
		//generate external utility
		writer.close(); // we close the file
	}
}

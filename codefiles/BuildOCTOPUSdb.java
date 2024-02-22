/*
del *.class
javac -cp ".;octopus_jars/*" BuildOCTOPUSdb.java
 java -Xmx32G -cp ".;octopus_jars/*" BuildOCTOPUSdb BVBRC_genomes.KmerSpecies.txt.gz -g 3
 java -Xmx32G -cp ".;octopus_jars/*" BuildOCTOPUSdb who_priority_bacteria.KmerSpecies.txt.gz -g 2
 java -Xmx32G -cp ".;octopus_jars/*" BuildOCTOPUSdb viral.1.1.genomic.KmerSpecies.txt.gz
 java -Xmx32G -cp ".;octopus_jars/*" BuildOCTOPUSdb mitochondrion.1.1.genomic.KmerSpecies.txt.gz
 java -Xmx32G -cp ".;octopus_jars/*" BuildOCTOPUSdb megares_database_v3.00.KmerSpecies.txt.gz
*/

import java.io.*;
import java.lang.*;
import java.math.*;
import java.nio.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.zip.GZIPInputStream;
import java.util.concurrent.TimeUnit;

import com.clearspring.analytics.hash.MurmurHash;
import com.greplin.bloomfilter.*;
import com.greplin.bloomfilter.allocator.*;
import org.h2.mvstore.*;
import org.mapdb.*;
import org.mapdb.volume.*;

public class BuildOCTOPUSdb
{
	public static long createChunk(HashMap<Integer,Long> sppKmer, BufferedReader r, LinkedList<Integer> seeds, BufferedWriter bw_tobehashed, String datafilerootname, boolean android, double skipprob) throws IOException
	{
		System.gc();
		float allram = (float)(Runtime.getRuntime().maxMemory());
		int seed = (int)seeds.removeFirst();
		int maxstoresize = 2147483646;
		double bffpr = 0.025d;
		TreeMap<Integer,Integer> ks_store = new TreeMap();
		boolean maxmemornum = false;
		FileWriter fw_justhashed = new FileWriter(datafilerootname+"/"+"justhashed.txt");
		BufferedWriter bw_justhashed = new BufferedWriter(fw_justhashed);
		String line = r.readLine();
		while(line!=null && !line.equals("") && !line.equals("\n") && !line.equals("\r\n"))
		{
			String [] datum = line.split(",");
			long kmer = Long.parseLong(datum[0]);
			String [] species_set = datum[1].split(";");
			if ( Math.random() < (double)(1d/(double)(species_set.length)) && Math.random() >= skipprob)
			//if (true)
			{
				int chosen = (int)(Math.random()*species_set.length);
				int spec = Integer.parseInt(species_set[chosen]);
				byte[] byteslongkmer = ByteBuffer.allocate(8).putLong(kmer).array();
				if (seeds.size()==0)
				{
					bw_tobehashed.write(kmer+","+spec); bw_tobehashed.newLine();
				}
				else
				{
					int hkmer = MurmurHash.hash(byteslongkmer,seed);
					if (ks_store.get(hkmer)==null)
					{
						ks_store.put(hkmer,spec);
						bw_justhashed.write(kmer+""); bw_justhashed.newLine();
						Long ski = sppKmer.get(spec);
						if (ski==null) {sppKmer.put(spec,1l);} else {sppKmer.put(spec,ski+1l);}
					}
					else
					{
						bw_tobehashed.write(kmer+","+spec); bw_tobehashed.newLine();
					}
				}
			}
			float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
			if (Math.random()<0.0000001) {System.out.println("\t\t"+ks_store.size()+" elements hashed ("+100*usedram/allram+"% heap memory used)");}
			if (usedram/allram>0.95)
			{
				System.out.println("\tMax resources (memory or integer hash segment) reached: ");
				System.out.println("\t\t"+100*usedram/allram+"% heap memory used");
				maxmemornum = true;
				break;	
			}
			line = r.readLine();
		}
		long size = (long)ks_store.size();
		if (size==0) {bw_justhashed.close(); fw_justhashed.close(); Files.deleteIfExists(Paths.get(datafilerootname+"/"+"justhashed.txt")); return 0;}
		System.out.println("\tSize of local store is "+size+" elements"); //System.out.println("seed "+seed);
		if (!android)
		{
			System.out.println("\t\tFlushing data on to file (Bloom filter and Mmap files)..");
			Volume volu = MappedFileVol.FACTORY.makeVolume(datafilerootname+"/int"+seed+".db",false);
			SortedTableMap.Sink<Integer,Integer> sink = SortedTableMap.create(volu,Serializer.INTEGER,Serializer.INTEGER).createFromSink();
			while(ks_store.size()>0) {Map.Entry<Integer,Integer> entry = ks_store.pollFirstEntry(); sink.put(entry.getKey(),entry.getValue());}
			sink.create(); volu.close();
		}
		else
		{
			System.out.println("\t\tFlushing data on to file (Bloom filter and Btree files)..");
			MVStore mvs = new MVStore.Builder().fileName(datafilerootname+"/int"+seed+".db").compress().cacheSize(128).pageSplitSize(64).open(); //16; 
			MVMap<Integer,Integer> bmap = mvs.openMap("mvs");
			while(ks_store.size()>0) {Map.Entry<Integer,Integer> entry = ks_store.pollFirstEntry(); bmap.put(entry.getKey(),entry.getValue());}
			mvs.close();
		}
		//BloomFilter<Long> bf = BloomFilter.create(Funnels.longFunnel(),(int)size,bffpr);
		File onDiskFile = new File(datafilerootname+"/"+"blo"+seed+".bl");
		BloomFilter bf = new BloomFilter.NewBuilder(onDiskFile,(int)size,bffpr).bucketSize(BucketSize.ONE).force(true).build();
		BufferedReader ri = new BufferedReader(new FileReader(datafilerootname+"/"+"justhashed.txt"),16384);
		String linei = ri.readLine();
		while(linei!=null && !linei.equals("") && !linei.equals("\n") && !linei.equals("\r\n"))
		{
			byte[] byteslongkmer = ByteBuffer.allocate(8).putLong(Long.parseLong(linei)).array();
			bf.add(byteslongkmer);
			linei = ri.readLine();
		}
		//printBloomFilterToFile(datafilerootname+"/"+"blo"+seeds[round]+".bl",bf);
		ri.close(); bf.close(); bw_justhashed.close(); fw_justhashed.close();
		Files.deleteIfExists(Paths.get(datafilerootname+"/"+"justhashed.txt"));
		ks_store = null;
		System.gc();
		if (maxmemornum) {System.out.println("\tMin. perfect hash chunk filled\r\n\t\tPassing on to next segment.."); size+=createChunk(sppKmer,r,seeds,bw_tobehashed,datafilerootname,android,skipprob);}
		return size;
	}
	
	public static void main(String[] args) throws IOException
	{
		long time0 = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;

		String datafile = "";
		boolean android = false;
		int maxgb = -1;
		double skipprob = 0;
		if (args==null || args.length==0 || (args.length>0 && args[0].startsWith("-h")))
		{
			System.out.println("Please run the program as"); 
			System.out.println("\t(Windows) java -Xmx[desired_ram] -cp \".;octopus_jars/*\" BuildOCTOPUSdb input_file [-a, -g GB]");
			System.out.println("\t(UNIX/Linux) java -Xmx[desired_ram] -cp \".:octopus_jars/*\" BuildOCTOPUSdb input_file [-a, -g GB]");
			System.out.println("\tThe input_file is the one generated by GenomesToKmers");
			System.out.println("\tThe -a option will create a database for the Android version");
			System.out.println("\tThe -g option will shrink the database to the desired GB size");
			System.out.println("\tThe desired_ram should be ~5G for datasets up to 50 million kmers, ~10G up to 100 million kmers, ~20G for up to 200 million kmers, etc. Even with suboptimal RAM allocation, it will work with any kmer cardinality but the minimal perfect hashing will be less efficient.");
			System.out.println("\tTyping java BuildOCTOPUSdb -h or -help will print this help");
			System.exit(0);
		}
		if (args!=null && args.length==1) datafile = args[0];
		if (args!=null && args.length>1)
		{
			for (int l=0; l<args.length; l++)
			{
				if (args[l].startsWith("-a")) {android=true;}
				else
				{
					if (args[l].startsWith("-g"))
					{
						l=l+1;
						maxgb=Integer.parseInt(args[l]);
						if (maxgb<1) System.out.println("Invalid GB size: will use all data");
					}
					else {if (!args[l].startsWith("-")) datafile = args[l];}
				}
			}
		}
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		System.out.println("Max RAM allocated (Xmx): "+allram/1048576+" MB");
		
		String datafilerootname = datafile.substring(Math.max(0,datafile.lastIndexOf("/")+1));
		datafilerootname = datafile.substring(Math.max(0,datafile.lastIndexOf("\\")+1));
		datafilerootname = datafilerootname.substring(0,Math.min(datafilerootname.indexOf(".KmerSpecies."),datafilerootname.length()));
		datafilerootname = datafilerootname + "_OCTOPUSdb";
		if (android) datafilerootname = datafilerootname + "_Android";
		System.out.println("Setting up directory "+datafilerootname);
		
		Integer[] seedA = {42, 69, 66, 77, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199};
		LinkedList<Integer> seeds = new LinkedList<Integer>(); Collections.addAll(seeds,seedA);
		
		System.out.println("Reading kmer and species file while");
		System.out.println("\tperforming minimal perfect hashing (from long to integer, max "+seedA.length+" rounds)");
		if (!android) {System.out.println("\t\tand storing hashes into memory-mapped files (MapDB)");}
		else {System.out.println("\t\tand storing hashes into b-tree files (MVStore)");}
		startTime = System.currentTimeMillis();
		BufferedReader r = new BufferedReader(new FileReader(datafile),DEFAULT_BUFFER_SIZE);
		File dbdir = new File(datafilerootname+"/"); if (!dbdir.exists()){dbdir.mkdirs();}
		for(File file: dbdir.listFiles()) if (!file.isDirectory()) file.delete();
		if (datafile.endsWith(".gz")) {r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(datafile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
		String line = r.readLine();
		if (!line.startsWith("species=")) {System.out.println("Wrong file format"); System.exit(0);}
		int s = Integer.parseInt(line.split("=")[1]);
		HashMap<Integer,String> spp_name = new HashMap(s);
		HashMap<Integer,Long> sppKmer = new HashMap(100000);
		while(true)
		{
			line = r.readLine();
			if (line==null || line.equals("") || line.equals("\n") || line.equals("\r\n")) {System.out.println("Wrong file format"); System.exit(0);}
			if (line.startsWith("k=")) {break;}
			int s_num = Integer.parseInt(line.split(",")[0]);
			String s_nam = line.split(",")[1];
			spp_name.put(s_num,s_nam);
		}
		if (line==null || line.equals("") || line.equals("\n") || line.equals("\r\n") || !line.startsWith("k=")) {System.out.println("Wrong file format"); System.exit(0);}
		int k = Integer.parseInt(line.split("=")[1]);
		line = r.readLine();
		if (line==null || line.equals("") || line.equals("\n") || line.equals("\r\n") || !line.startsWith("mask_nonACGTchars=")) {System.out.println("Wrong file format"); System.exit(0);}
		boolean mask = Boolean.parseBoolean(line.split("=")[1]);
		line = r.readLine();
		if (line==null || line.equals("") || line.equals("\n") || line.equals("\r\n") || !line.startsWith("number_kmers=")) {System.out.println("Wrong file format"); System.exit(0);}
		long n1 = Long.parseLong(line.split("=")[1]);
		System.out.println("Initial number of kmers: "+n1+" (including kmers shared by multiple species)");
		if (maxgb>0) skipprob = 1d-(double)maxgb/(n1*0.00000000802d);
		long nfinal = 0;
		while(true)
		{
			FileWriter fw_tobehashed = new FileWriter(datafilerootname+"/"+"tobehashed.txt");
			BufferedWriter bw_tobehashed = new BufferedWriter(fw_tobehashed);
			nfinal += createChunk(sppKmer,r,seeds,bw_tobehashed,datafilerootname,android,skipprob);
			bw_tobehashed.close(); fw_tobehashed.close(); r.close();
			Path source = Paths.get(datafilerootname+"/"+"tobehashed.txt"); 
			if (Files.size(source)==0 || seeds.size()==0) {break;}
			Files.move(source,source.resolveSibling("source.txt"),StandardCopyOption.REPLACE_EXISTING);
			r = new BufferedReader(new FileReader(datafilerootname+"/"+"source.txt"),DEFAULT_BUFFER_SIZE);
			System.gc();
			System.out.println("\tSource file completed, checking remaining kmers to be hashed ("+Files.size(Paths.get(datafilerootname+"/"+"source.txt"))/(1024*1024)+" MB to go)");
		}
		r.close();
		Files.deleteIfExists(Paths.get(datafilerootname+"/"+"source.txt"));
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("Time employed to finalize minimal perfect hashing and memory mapped files: "+elapsedTime/1000+" seconds");
		
		boolean longdb = false;
		File file = new File(datafilerootname+"/"+"tobehashed.txt");
		if (file.length()!=0)
		{
			System.out.print("Checking for leftover (long) hashes..");
			longdb = true;
			r = new BufferedReader(new FileReader(datafilerootname+"/"+"tobehashed.txt"),DEFAULT_BUFFER_SIZE);
			line = r.readLine();
			if (!android)
			{
				Volume volu = MappedFileVol.FACTORY.makeVolume(datafilerootname+"/long.db",false);
				SortedTableMap.Sink<Long,Integer> sink = SortedTableMap.create(volu,Serializer.LONG,Serializer.INTEGER).createFromSink();
				while(line!=null)
				{
					if (!line.equals("") && !line.equals("\n") && !line.equals("\r\n"))
					{
						String [] datum = line.split(",");
						long kmer = Long.parseLong(datum[0]);
						int spec = Integer.parseInt(datum[1]);
						sink.put(kmer,spec);
						nfinal++;
					}
					line = r.readLine();
				}
				sink.create(); volu.close();
			}
			else
			{
				MVStore mvs = new MVStore.Builder().fileName(datafilerootname+"/long.db").compress().cacheSize(128).pageSplitSize(64).open(); //16; 
				MVMap<Long,Integer> bmap = mvs.openMap("mvs");
				while(line!=null)
				{
					if (!line.equals("") && !line.equals("\n") && !line.equals("\r\n"))
					{
						String [] datum = line.split(",");
						long kmer = Long.parseLong(datum[0]);
						int spec = Integer.parseInt(datum[1]);
						bmap.put(kmer,spec);
						nfinal++;
					}
					line = r.readLine();
				}
				mvs.close();
			}
			r.close();
			System.out.println("done.");
		}
		Files.deleteIfExists(Paths.get(datafilerootname+"/"+"tobehashed.txt"));
		System.out.println("Final number of kmers: "+nfinal+" (after probabilistic filtering based on species' multiplicity)");
		
		File directoryPath = new File(datafilerootname);
		String fileList[] = directoryPath.list();
		seeds = new LinkedList();
		for(int i=0; i<fileList.length; i++) {if (fileList[i].indexOf("blo")!=-1) {seeds.add(Integer.parseInt(fileList[i].replaceAll("[^0-9]","")));}}
		
		System.out.println("Calculating classification thresholds");
		Volume[] sm;
		SortedTableMap<Integer,Integer>[] mm;
		sm = new Volume[seeds.size()];
		mm = new SortedTableMap[seeds.size()];
		Volume longvm = null;
		SortedTableMap<Long,Integer> longmapm = null;
		MVStore[] sa;
		MVMap<Integer,Integer>[] ma;
		sa = new MVStore[seeds.size()];
		ma = new MVMap[seeds.size()];
		MVStore longva = null;
		MVMap<Integer,Integer> longmapa = null;
		BloomFilter[] bf;
		if (!android)
		{
			for (int i=0; i<seeds.size(); i++)
			{
				sm[i] = MappedFileVol.FACTORY.makeVolume(datafilerootname+"/int"+seeds.get(i)+".db",true);
				mm[i] = SortedTableMap.open(sm[i],Serializer.INTEGER,Serializer.INTEGER);
			}
			bf = new BloomFilter[seeds.size()];
			for (int i=0; i<seeds.size(); i++) {bf[i] = new BloomFilter.OpenBuilder(new File(datafilerootname+"/blo"+seeds.get(i)+".bl")).build();}
			if (longdb)
			{
				longvm = MappedFileVol.FACTORY.makeVolume(datafilerootname+"/long.db",true);
				longmapm = SortedTableMap.open(longvm,Serializer.LONG,Serializer.INTEGER);
			}
		}
		else
		{
			for (int i=0; i<seeds.size(); i++)
			{
				sa[i] = new MVStore.Builder().fileName(datafilerootname+"/int"+seeds.get(i)+".db").readOnly().open();
				ma[i] = sa[i].openMap("mvs");
			}
			bf = new BloomFilter[seeds.size()];
			for (int i=0; i<seeds.size(); i++) {bf[i] = new BloomFilter.OpenBuilder(new File(datafilerootname+"/blo"+seeds.get(i)+".bl")).build();}
			longmapa = null;
			longva = null;
			if (longdb)
			{
				longva = new MVStore.Builder().fileName(datafilerootname+"/long.db").readOnly().open();
				longmapa = longva.openMap("mvs");
			}
		}
		int[] readlengths = {1,32,38,45,54,64,76,91,108,128,152,181,215,256,304,362,431,512,609,724,861,1024,1218,1448,1722,2048,2435,2896,3444,4096,4871,5793,6889,8192,9742,11585,13777,16384,19484,23170,27554,32768,38968,46341,55109,65536,77936,92682,110218,131072,155872,185364,220436,262144,311744,370728,440872,524288,623487,741455,881744,1048576,1246974,1482910,1763488,2097152,2493948,2965821,3526975,4194304,4987896,5931642,7053950,8388608};
		System.out.println("\tRead lengths range: "+readlengths[0]+" to "+readlengths[readlengths.length-1]+" nucleotide bases");
		readlengths[0] = k-4;
		double er = 0d; if (mask) er = 0.00005d; 
		int alphabetSizePrime = 5;
		if (!mask) alphabetSizePrime = 4;
		long [][] ppn = precomputedPowsForNucs(alphabetSizePrime,33);
		HashMap<Character,Integer> nc = nucleotideCodes(mask);
		HashMap<Character,Integer> ncr = nucleotideCodesRev(mask);
		int boot = 1000;
		
		int step = readlengths[0];
		float estmem = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		for (int j=1; j<readlengths.length; j++)
		{
			step = readlengths[j]-readlengths[j-1];
			char[] cazzo = randomSequence(step+k-4,er);
			TreeSet<Long> th = hashSequence(cazzo, k-4, ppn, nc, ncr, mask);
			estmem = Math.max(estmem,(float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
			if (estmem/allram>0.9123f) {break;}
		}
		double bffpr = 0.0001d;
		HashMap<Integer,ArrayList<Integer>> bootreadhits = new HashMap(readlengths.length);
		for (int j=0; j<readlengths.length; j++) {bootreadhits.put(readlengths[j],new ArrayList());}
		int readlengthslimit = readlengths.length;
		System.out.print("\t\t");
		for (int b=0; b<boot; b++)
		{
			if (b>4) {readlengthslimit = Arrays.binarySearch(readlengths,4194304);}
			if (b>8) {readlengthslimit = Arrays.binarySearch(readlengths,2097152);}
			if (b>16) {readlengthslimit = Arrays.binarySearch(readlengths,1048576);}
			if (b>32) {readlengthslimit = Arrays.binarySearch(readlengths,524288);}
			if (b>64) {readlengthslimit = Arrays.binarySearch(readlengths,262144);}
			if (b>128) {readlengthslimit = Arrays.binarySearch(readlengths,131072);}
			if (b>256) {readlengthslimit = Arrays.binarySearch(readlengths,65536);}
			if (b>512) {readlengthslimit = Arrays.binarySearch(readlengths,32768);}
			if (readlengthslimit<0) {readlengthslimit = readlengths.length;}
			BloomFilter readhashes = new BloomFilter.NewBuilder(null,readlengths[readlengths.length-1],bffpr).bucketSize(BucketSize.ONE).force(true).build();
			TreeMap<Integer,Integer> sppHits = new TreeMap();
			char[] kmerremainder = randomSequence(k-4,er);
			TreeSet<Long> readhashestemp = hashSequence(kmerremainder, k-4, ppn, nc, ncr, mask);
			byte[] readba = ByteBuffer.allocate(8).putLong(readhashestemp.pollFirst()).array();
			readhashes.add(readba);
			for (int p=0; p<readlengthslimit; p++)
			{
				int counthits = 0;
				int readlengthtemp = 0;
				int readlengthtotal = k-4;
				if (p>0) {readlengthtotal = readlengthtotal + readlengths[p] - readlengths[p-1];}
				while (readlengthtemp<readlengthtotal)
				{
					char[] read = randomSequence(Math.min(readlengthtotal,Math.min(readlengthtotal-readlengthtemp,step)),er);
					for (int j=0; j<kmerremainder.length; j++) {read[j]=kmerremainder[j];}
					for (int j=0; j<kmerremainder.length; j++) {kmerremainder[j]=read[read.length-k+4+j];}
					readhashestemp = hashSequence(read, k-4, ppn, nc, ncr, mask);
					while (readhashestemp.size()>0)
					{
						long ah = readhashestemp.pollFirst();
						byte[] lb = ByteBuffer.allocate(8).putLong(ah).array();
						if (!readhashes.contains(lb))
						{
							readhashes.add(lb);
							Integer spp = null;
							for (int i=0; i<seeds.size(); i++)
							{
								if (bf[i].contains(lb))
								{
									int mh = MurmurHash.hash(lb,seeds.get(i));
									if (!android) {spp = mm[i].get(mh);} else {spp = ma[i].get(mh);} 
									if (spp!=null) {break;}
								}
							}
							if (longdb && spp==null)
							{
								if (!android) {spp = longmapm.get(ah);} else {spp = longmapa.get(ah);}
							}
							if (spp!=null)
							{
								Integer ch = sppHits.get(spp);
								if (ch!=null) {sppHits.put(spp,ch+1);} else {sppHits.put(spp,1);}
							}
						}
					}
					for (Map.Entry<Integer,Integer> entry : sppHits.entrySet()) {counthits = (int)(Math.max(counthits,entry.getValue()));}
					counthits = (int)Math.round(counthits+counthits*bffpr);
					readlengthtemp = readlengthtemp + read.length;
				}
				ArrayList<Integer> ar = bootreadhits.get(readlengths[p]); ar.add(counthits); bootreadhits.put(readlengths[p],ar);
			}
			readhashes.close();
			if (boot>10) {if (b%(boot/10)==0) {System.out.print((int)(100d*((double)b/boot))+"%.. ");}}
		}
		System.out.println("done");
		ArrayList<Double> hits = new ArrayList();
		ArrayList<Double> hits_variance = new ArrayList();
		for (int j=0; j<readlengths.length; j++)
		{
			ArrayList<Integer> hitstemp = bootreadhits.get(readlengths[j]);
			//Collections.sort(hitstemp);
			//int expected = hitstemp.get((int)(hitstemp.size()*0.50));
			double expected=0; for (int p=0; p<hitstemp.size(); p++) {expected+=hitstemp.get(p);} expected=expected/hitstemp.size();
			hits.add(expected);
			double variance=0; for (int p=0; p<hitstemp.size(); p++) {variance=(expected-hitstemp.get(p))*(expected-hitstemp.get(p));} variance=variance/hitstemp.size();
			hits_variance.add(variance);
		}
		Collections.sort(hits);
		Collections.sort(hits_variance);
		System.out.println("Expected average (variance) random hits per species for read lengths of:");
		System.out.println("\t"+readlengths[0]+" = "+Math.round(hits.get(0)*10000)/10000d+" ("+Math.round(hits_variance.get(0)*10000)/10000d+") from "+bootreadhits.get(readlengths[0]).size()+" bootstraps");
		System.out.println("\t"+readlengths[readlengths.length/4]+" = "+Math.round(hits.get(readlengths.length/4)*10000)/10000d+" ("+Math.round(hits_variance.get(readlengths.length/4)*10000)/10000d+") from "+bootreadhits.get(readlengths[readlengths.length/4]).size()+" bootstraps");
		System.out.println("\t"+readlengths[readlengths.length/2]+" = "+Math.round(hits.get(readlengths.length/2)*10000)/10000d+" ("+Math.round(hits_variance.get(readlengths.length/2)*10000)/10000d+") from "+bootreadhits.get(readlengths[readlengths.length/2]).size()+" bootstraps");
		System.out.println("\t"+readlengths[3*readlengths.length/4]+" = "+Math.round(hits.get(3*readlengths.length/4)*10000)/10000d+" ("+Math.round(hits_variance.get(3*readlengths.length/4)*10000)/10000d+") from "+bootreadhits.get(readlengths[3*readlengths.length/4]).size()+" bootstraps");
		System.out.println("\t"+readlengths[readlengths.length-1]+" = "+Math.round(hits.get(readlengths.length-1)*10000)/10000d+" ("+Math.round(hits_variance.get(readlengths.length-1)*10000)/10000d+") from "+bootreadhits.get(readlengths[readlengths.length-1]).size()+" bootstraps");
		
		if (!android) {for (int i=0; i<sm.length; i++) sm[i].close(); if (longvm!=null) longvm.close();} else {for (int i=0; i<sa.length; i++) sa[i].close(); if (longva!=null) longva.close();}
		for (int i=0; i<seeds.size(); i++) bf[i].close();
		
		System.out.println("Writing metadata information on to file");
		System.out.print("\t");
		FileWriter fw = new FileWriter(datafilerootname+"/"+"info.txt");
		BufferedWriter bw = new BufferedWriter(fw);
		bw.write("k,"+k); bw.newLine();
		bw.write("mask_nonACGTchars,"+mask); bw.newLine();
		bw.write("number_kmers,"+nfinal); bw.newLine();
		bw.write("randomhits_thresholds,");
		for (int p=0; p<hits.size()-1; p++) {bw.write(readlengths[p]+":"+hits.get(p)+";");}
		bw.write(readlengths[hits.size()-1]+":"+hits.get(hits.size()-1)); bw.newLine();
		bw.write("long_key_value_store,"+longdb); bw.newLine();
		bw.write("int_key_value_stores_seeds,"); for(int q=0; q<seeds.size()-1; q++) {bw.write(seeds.get(q)+";");} bw.write(seeds.get(seeds.size()-1)+""); bw.newLine();
		Collection<Integer> keysc = sppKmer.keySet();
		ArrayList<Integer> keys = new ArrayList<Integer>(keysc);
		int gg=0;
		for (int key : keys)
		{
			long counts = sppKmer.get(key);
			int si = key;
			long ki = counts;
			String ni = spp_name.get(si);
			bw.write(si+","+ki+","+ni); bw.newLine();
			if (keys.size()>=15)
			{
				if ((gg%(keys.size()/6)==0))
				{
					System.out.print(100*gg/keys.size()+"%.. ");
				}
			}
			gg++;
		}
		bw.flush(); bw.close(); fw.close(); r.close();
		System.out.println("done.");
		
		elapsedTime = System.currentTimeMillis() - time0;
		System.out.println("Total time employed "+elapsedTime/1000+" seconds");
	}
	
	public static char[] randomSequence(int n, double error)
	{
		char[] k = new char[n];
		for (int i=0; i<n; i++)
		{
			double d = Math.random();
			if (d<0.25d) {k[i]='A';}
				else if (d<0.50d) {k[i]='C';}
					else if (d<0.75d) {k[i]='G';}
						else {k[i]='T';}
			if (Math.random()<error) {k[i]='N';}
		}
		return k;
	}
	
	public static HashMap<Character,Integer> nucleotideCodes(boolean mask)
	{
		HashMap<Character,Integer> hs = new HashMap(23);
		if ( mask) {hs.put('A',0); hs.put('a',0); hs.put('C',1); hs.put('c',1); hs.put('G',2); hs.put('g',2); hs.put('N',3); hs.put('n',3); hs.put('T',4); hs.put('t',4); hs.put('U',4); hs.put('u',4); }
		if (!mask) {hs.put('A',0); hs.put('a',0); hs.put('C',1); hs.put('c',1); hs.put('G',2); hs.put('g',2); hs.put('T',3); hs.put('t',3); hs.put('U',3); hs.put('u',3); }
		return hs;
	}
	public static HashMap<Character,Integer> nucleotideCodesRev(boolean mask)
	{
		HashMap<Character,Integer> hs = new HashMap(23);
		if ( mask) {hs.put('A',4); hs.put('a',4); hs.put('C',2); hs.put('c',2); hs.put('G',1); hs.put('g',1); hs.put('N',3); hs.put('n',3); hs.put('T',0); hs.put('t',0); hs.put('U',0); hs.put('u',0); }
		if (!mask) {hs.put('A',3); hs.put('a',3); hs.put('C',2); hs.put('c',2); hs.put('G',1); hs.put('g',1); hs.put('T',0); hs.put('t',0); hs.put('U',0); hs.put('u',0); }
		return hs;
	}
	public static long [][] precomputedPowsForNucs(int p, int k)
	{
		long [][] ppn = new long [p][];
		for (int w=0; w<p; w++)
		{
			long [] pp = new long [k];
			for (int i=0; i<k; i++)
			{
				long pow = 1;
				for (int j=0; j<i; j++)
					pow = pow*p;
				pp[i] = pow;
			}
			for (int i=0; i<k; i++) pp[i]=pp[i]*w;
			ppn[w]=pp;
		}
		return ppn;
	}
	public static long polynomialHash(int[] str, int start, int stop, long [][] pp)
	{
		long hash_val = 0;
		int k = stop - start;
		for (int i=0; i<k; i++)
		{
			hash_val = hash_val + pp[str[(stop-1)-i]][i];
		}
		return hash_val;
	}
	public static long nextHash(long [][] pp, long prev_hash, int k, int first, int next)
	{
		//long nextH = ( prev_hash - pp[first][k-1] ) * pp[1][1] + next;
		long nextH = ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + next;
		if (pp.length==5) nextH += prev_hash - pp[first][k-1];
		return nextH;
	}
	public static long getMinimizer(long[] kmers, long mink, int start, int step)
	{
		long minimizer = Long.MAX_VALUE;
		if (start==0)
		{
			boolean found = false;
			for (int i=0; i<=step; i++) {if (kmers[i]>=0 && kmers[i]<=minimizer) {found = true; minimizer = kmers[i];}}
			if (found) {return minimizer;} else {return -1l;}
		}
		else
		{
			long oldk = kmers[start-1];
			long newk = kmers[start+step];
			if (mink==-1) return newk;
			if (oldk==-1) return (long)Math.min(mink,newk); //mink cannot be -1 here
			if (mink==oldk)  //neither mink nor oldk can be -1 here
			{
				boolean found = false;
				for (int i=start; i<=(start+step); i++) {if (kmers[i]>=0 && kmers[i]<=minimizer) {found = true; minimizer = kmers[i];}}
				if (found) {return minimizer;} else {return -1l;}
			}
			if (newk>=0) {return (long)Math.min(mink,newk);} else {return mink;} 
		}
	}
	public static TreeSet<Long> hashSequence(char[] fwdCharArray, int k, long [][] pp, HashMap<Character,Integer> nc, HashMap<Character,Integer> ncr, boolean mask)
	{
		int[] fwdIntArray = new int[fwdCharArray.length];
		int[] rwdIntArray = new int[fwdCharArray.length];
		TreeSet<Integer> badFwd = new TreeSet();
		for (int g=0; g<fwdCharArray.length; g++) {Integer ic = nc.get(fwdCharArray[g]); if (ic!=null) {fwdIntArray[g] = ic;} else {if (mask) {fwdIntArray[g] = 3;} else {fwdIntArray[g] = 0; ; for (int q=-(k-1); q<=0; q++) {badFwd.add(g+q);}}}}
		for (int g=0; g<fwdCharArray.length; g++) {Integer ic = ncr.get(fwdCharArray[fwdCharArray.length-(g+1)]); if (ic!=null) {rwdIntArray[g] = ic;} else {if (mask) {rwdIntArray[g] = 3;} else {rwdIntArray[g] = 3;}}}
		long fh=polynomialHash(fwdIntArray,0,k,pp);
		long rh=polynomialHash(rwdIntArray,0,k,pp);
		long [] fwd_hashes = new long [fwdCharArray.length-k+1]; fwd_hashes[0]=fh;
		long [] rwd_hashes = new long [fwdCharArray.length-k+1]; rwd_hashes[fwdCharArray.length-k]=rh;
		for (int g=1; g<fwdCharArray.length-k+1; g++)
		{
			int first = fwdIntArray[g-1];
			int next = fwdIntArray[g+k-1];
			fh = nextHash(pp,fh,k,first,next);
			fwd_hashes[g]=fh;
			first = rwdIntArray[g-1];
			next = rwdIntArray[g+k-1];
			rh = nextHash(pp,rh,k,first,next);
			rwd_hashes[fwdCharArray.length-k-g]=rh;
		}
		for (int g=0; g<fwd_hashes.length; g++)
		{
			if (!badFwd.contains(g)) {if (rwd_hashes[g]<fwd_hashes[g]) fwd_hashes[g]=rwd_hashes[g];}
			else {fwd_hashes[g]=-1;}
		}
		int step = 4;
		//System.out.println(fwd+"\r\n"+"original hashes = "+Arrays.toString(fwd_hashes));
		TreeSet<Long> minimizer_hashes = new TreeSet();
		if ((fwd_hashes.length-step)>0)
		{
			long min = Long.MAX_VALUE;
			for (int g=0; g<fwd_hashes.length-step; g++)
			{
				long minT = getMinimizer(fwd_hashes,min,g,step);
				minimizer_hashes.add(minT);
				min = minT;
			}
		}
		else
		{
			boolean found = false;
			long min = Long.MAX_VALUE;
			for (int g=0; g<fwd_hashes.length; g++)
			{
				if (fwd_hashes[g]>0 && fwd_hashes[g]<min) {min=fwd_hashes[g]; found=true;}
			}
			if (found) {minimizer_hashes.add(min);} else  {minimizer_hashes.add(-1l);}
		}
		//System.out.println("minimize hashes = "+Arrays.toString(minimizer_hashes));
		return minimizer_hashes;
	}

}

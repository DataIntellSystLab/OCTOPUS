//Oligomer-based Classifier of Taxonomic Operational and Pan-genome Units via Singletons (OCTOPUS)

/*
del *.class
javac -cp ".;octopus_jars/*" OCTOPUS.java
 java -cp ".;octopus_jars/*" OCTOPUS d:who_priority_bacteria_OCTOPUSdb f:..\simulmix.fastq
 java -cp ".;octopus_jars/*" OCTOPUS d:viral.1.1.genomic_OCTOPUSdb f:..\simulmix.fastq t:16 s:0.75
 java -cp ".;octopus_jars/*" OCTOPUS d:mitochondrion.1.1.genomic_OCTOPUSdb f:..\sra_data.fastq.gz t:16
 java -cp ".;octopus_jars/*" OCTOPUS d:megares_database_v3.00_OCTOPUSdb f:..\simulmix.fastq t:16
*/

import java.io.*;
import java.io.File;
import java.lang.*;
import java.math.*;
import java.nio.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.concurrent.*;
import java.util.concurrent.atomic.*;
import org.mapdb.*;
import org.mapdb.volume.*;
import com.clearspring.analytics.stream.cardinality.*;
import com.clearspring.analytics.hash.MurmurHash;
import com.greplin.bloomfilter.*;
import com.greplin.bloomfilter.allocator.*;

class Tentacle extends Thread
{
	public int nameT;
	public ConcurrentLinkedQueue<ClassificationResult> resultsT;
	public String headerT;
	public String fwdT;
	public int kT;
	public boolean maskT;
	public long[][] ppnT;
	public HashMap<Character,Integer> ncT;
	public HashMap<Character,Integer> ncrT;
	public KmerMetaDataBase storeT;
	
	public Tentacle(int name, ConcurrentLinkedQueue<ClassificationResult> results, String header, String fwd, int k, boolean mask, long [][] ppn, HashMap<Character,Integer> nc, HashMap<Character,Integer> ncr, KmerMetaDataBase store) 
	{
		nameT=name;
		resultsT=results;
		headerT=header;
		fwdT=fwd;
		kT=k;
		maskT=mask;
		ppnT=ppn;
		ncT=nc;
		ncrT=ncr;
		storeT=store;
	}
   
	public void run()
	{
		try
		{
			ClassificationResult cr = new ClassificationResult();
			cr.read_header = headerT;
			HashMap<Long,Integer> hits_posT = new HashMap();
			TreeSet<Long> fwd_hashes = new TreeSet();
			int read_kmers = 0;
			if (fwdT.length()>=kT)
			{
				fwd_hashes = OCTOPUS.hashSequence(fwdT, kT, ppnT, ncT, ncrT, maskT);
				read_kmers = fwd_hashes.size();
				while(fwd_hashes.size()>0)
				{
					long ah = fwd_hashes.pollFirst();
					if (ah>=0)
					{
						Integer h = storeT.lookUpKmerSpecies(ah); 
						if (h!=null) {hits_posT.put(ah,h);}
					}
				}
			}
			if (hits_posT.size()>0)
			{
				Collection<Long> keysc = hits_posT.keySet();
				ArrayList<Long> keys = new ArrayList<Long>(keysc);
				HashMap<Integer,Integer> sppf = new HashMap();
				for (long key : keys)
				{
					int spp = hits_posT.get(key);
					Integer sf = sppf.get(spp);
					if (sf!=null) {sppf.put(spp,sf+1);}
					else {sppf.put(spp,1);}
				}
				int smaxid=-1;
				int smaxfr=-1;
				Collection<Integer> keyscs = sppf.keySet();
				ArrayList<Integer> keyss = new ArrayList<Integer>(keyscs);
				for (int key : keyss)
				{
					int sf = sppf.get(key);
					if (sf>smaxfr) {smaxfr=sf;smaxid=key;}
				}
				Integer fk = storeT.ht.floorKey(fwdT.length());
				Integer ck = storeT.ht.ceilingKey(fwdT.length());
				int closestKey = fk;
				if (ck!=null) {if ( Math.abs(ck-fwdT.length()) < Math.abs(fk-fwdT.length()) ) closestKey = ck;}
				int thr = storeT.ht.get(closestKey);
				if ( smaxfr>thr )
				{
					LinkedList<Long> skh = new LinkedList();
					for (long key : keys)
					{
						int spp = hits_posT.get(key);
						if (spp==smaxid) {skh.add(key);}
					}
					cr = new ClassificationResult(headerT, read_kmers, smaxid, skh);
				}
			}
			resultsT.add(cr);
		}
		catch (Exception e) {System.out.println(e); System.exit(0);}
	}
}

class Species implements Serializable
{
	public long nkmers;
	public HyperLogLogPlus ukmersfound;
	public long readdepth;
	public Species()
	{
		nkmers=-1;
		ukmersfound=null;
		readdepth=-1;
	}
}

class ClassificationResult implements Serializable
{
	public String read_header;
	public int read_kmers;
	public int species_index;
	public LinkedList<Long> species_kmer_hits;
	public ClassificationResult()
	{
		read_header=null;
		read_kmers=0;
		species_index=-1;
		species_kmer_hits=null;
	}
	public ClassificationResult(String he, int rk, int si, LinkedList<Long> skh)
	{
		read_header=he;
		read_kmers=rk;
		species_index=si;
		species_kmer_hits=skh;
	}
}

class KmerMetaDataBase
{
	Volume[] s;
	SortedTableMap<Integer,Integer>[] m;
	BloomFilter[] bf;
	public boolean longdb;
	Volume longv;
	SortedTableMap<Long,Integer> longmap;
	public int k;
	public boolean mask;
	public long n;
	public TreeMap<Integer,Long> sk;
	public TreeMap<Integer,String> sn;
	public static int[] seeds;
	public TreeMap<Integer,Integer> ht;
	
	public KmerMetaDataBase()
	{
		s = null;
		m = null;
		bf = null;
		longdb = false;
		longv = null;
		longmap = null;
		k = 0;
		mask = false;
		n = 0;
		sk = null;
		sn = null;
		seeds = null;
		ht = null;
	}
	public KmerMetaDataBase(String ukmdb, double pvalue_or_minfreq) throws Exception
	{
		FileReader fr = new FileReader(ukmdb+"/info.txt");
		BufferedReader br = new BufferedReader(fr);
		k = Integer.parseInt(br.readLine().split(",")[1])-4;
		System.out.println("\tLength of kmers: "+(k+4)+"; length of minimizers: "+k);
		mask = Boolean.parseBoolean(br.readLine().split(",")[1]);
		System.out.println("\tNon-ACGT character masking with N set to: "+mask);
		n = Long.parseLong(br.readLine().split(",")[1]);
		ht = new TreeMap();
		String[] hits_thresholds = br.readLine().split(",")[1].split(";");
		if (pvalue_or_minfreq<1)
		{
			double stdev = invPhiCDF(pvalue_or_minfreq,-8d,8d);
			for (int g=0; g<hits_thresholds.length; g++) {String [] rlth = hits_thresholds[g].split(":"); double newthr = Double.parseDouble(rlth[1]); newthr = newthr+stdev*Math.sqrt(newthr); ht.put(Integer.parseInt(rlth[0]),(int)Math.round(newthr));}
		}
		else
		{ht.put(k,(int)(pvalue_or_minfreq-1)); ht.put(Integer.MAX_VALUE,(int)(pvalue_or_minfreq-1));}
		longdb = false;
		if (br.readLine().split(",")[1].equals("true")) {longdb = true;}
		String[] seedsvalues = br.readLine().split(",")[1].split(";");
		seeds = new int[seedsvalues.length];
		for (int i=0; i<seedsvalues.length; i++) {seeds[i] = Integer.parseInt(seedsvalues[i]);}
		sk = new TreeMap();
		sn = new TreeMap();
		String spnk = br.readLine();
		while (spnk!=null)
		{
			sk.put(Integer.parseInt(spnk.split(",")[0]),Long.parseLong(spnk.split(",")[1]));
			sn.put(Integer.parseInt(spnk.split(",")[0]),spnk.split(",")[2]);
			spnk = br.readLine();
		}
		br.close(); fr.close();
		s = new Volume[seeds.length];
		m = new SortedTableMap[seeds.length];
		for (int i=0; i<seeds.length; i++)
		{
			s[i] = MappedFileVol.FACTORY.makeVolume(ukmdb+"/int"+seeds[i]+".db",true);
			m[i] = SortedTableMap.open(s[i],Serializer.INTEGER,Serializer.INTEGER);
		}
		bf = new BloomFilter[seeds.length];
		for (int i=0; i<seeds.length; i++)
		{
			bf[i] = new BloomFilter.OpenBuilder(new File(ukmdb+"/blo"+seeds[i]+".bl")).build();
			//bf[i] = BloomFilter.readFrom(new FileInputStream(new File(ukmdb+"/blo"+seeds[i]+".bl")),Funnels.longFunnel());
			/*InputStream in1 = new FileInputStream(new File(ukmdb+"/blo"+seeds[i]+".bl"));
			ObjectInputStream oin1 = new ObjectInputStream(in1);
			bf[i] = (BloomFilter)oin1.readObject();
			in1.close();
			oin1.close();*/
		}
		longmap = null;
		longv = null;
		if (longdb)
		{
			longv = MappedFileVol.FACTORY.makeVolume(ukmdb+"/long.db",true);
			longmap = SortedTableMap.open(longv,Serializer.LONG,Serializer.INTEGER);
		}
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
	}
	public static double invPhiCDF(double prob, double lowerl, double higherl)
	{
		double dtol=0.00000001d;
		double midp=0.5*(higherl-lowerl)+lowerl;
		if (higherl-lowerl<dtol) {return midp;}
		double cdf=0;
		if (midp>8.0d)
		{
			cdf=1d;
		}
		else 
		{
			if(midp<-8.0d) 
			{
				cdf=0d;
			}
			else
			{ 
				double temp=midp;
				double midp2=Math.pow(midp,2);
				double cums=0d;
        			for (int j=3; temp+cums!=cums; j=j+2)
				{
           	 			cums+=temp;
           	 			temp=midp2*temp/j;
				}
				double pdf=Math.exp(-0.5d*midp2)/Math.sqrt(2*Math.PI);
				cdf=cums*pdf+0.5d;
			}
		}
		if (cdf>prob) {return invPhiCDF(prob,lowerl,midp);} else {return invPhiCDF(prob,midp,higherl);}
	}
	public int getKValue() {return k;}
	public boolean getMaskValue() {return mask;}
	public long getNumKmers() {return n;}
	public long getNumKmersSpecies(int key) {return sk.get(key);}
	public String getSpeciesName(int key) {return sn.get(key);}
	public void closeKmerMetaDataBase() throws IOException {for (int i=0; i<seeds.length; i++) {s[i].close(); bf[i].close();}; if(longdb){longv.close();}}
	public Collection<Integer> getSpeciesKmerKeySet() {return sk.keySet();}
	public Integer lookUpKmerSpecies(long ah)
	{
		Integer spp = null;
		byte[] lb = ByteBuffer.allocate(8).putLong(ah).array();
		for (int i=0; i<seeds.length; i++)
		{
			if (bf[i].contains(lb))
			{
				int mh = MurmurHash.hash(lb,seeds[i]);
				spp = m[i].get(mh);
				if (spp!=null) {return spp;}
			}
		}
		if (longdb) {spp = longmap.get(ah); if (spp!=null) {return spp;}}
		return spp;
	}
}

public class OCTOPUS
{
	public static long fastModulo(long dividend, long divisor)
	{
		if (((divisor-1) & divisor) == 0) {return (dividend & (divisor - 1));}
		return (dividend%divisor);
	}
	
	public static double arcsinhyp(double x)
	{
		final double s;
		if (Double.doubleToRawLongBits(x)<0)
		{
			x = Math.abs(x);
			s = -1.0d;
		}
		else
		{
			s = 1.0d;
		}
		return s * Math.log(Math.sqrt(x * x + 1.0d) + x);
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
	public static TreeSet<Long> hashSequence(String fwd, int k, long [][] pp, HashMap<Character,Integer> nc, HashMap<Character,Integer> ncr, boolean mask)
	{
		char[] fwdCharArray = fwd.toCharArray();
		int[] fwdIntArray = new int[fwdCharArray.length];
		int[] rwdIntArray = new int[fwdCharArray.length];
		TreeSet<Integer> badFwd = new TreeSet();
		for (int g=0; g<fwdCharArray.length; g++) {Integer ic = nc.get(fwdCharArray[g]); if (ic!=null) {fwdIntArray[g] = ic;} else {if (mask) {fwdIntArray[g] = 3;} else {fwdIntArray[g] = 0; ; for (int q=-(k-1); q<=0; q++) {badFwd.add(g+q);}}}}
		for (int g=0; g<fwdCharArray.length; g++) {Integer ic = ncr.get(fwdCharArray[fwdCharArray.length-(g+1)]); if (ic!=null) {rwdIntArray[g] = ic;} else {if (mask) {rwdIntArray[g] = 3;} else {rwdIntArray[g] = 3;}}}
		long fh=polynomialHash(fwdIntArray,0,k,pp);
		long rh=polynomialHash(rwdIntArray,0,k,pp);
		long [] fwd_hashes = new long [fwd.length()-k+1]; fwd_hashes[0]=fh;
		long [] rwd_hashes = new long [fwd.length()-k+1]; rwd_hashes[fwd.length()-k]=rh;
		for (int g=1; g<fwd.length()-k+1; g++)
		{
			int first = fwdIntArray[g-1];
			int next = fwdIntArray[g+k-1];
			fh = nextHash(pp,fh,k,first,next);
			fwd_hashes[g]=fh;
			first = rwdIntArray[g-1];
			next = rwdIntArray[g+k-1];
			rh = nextHash(pp,rh,k,first,next);
			rwd_hashes[fwd.length()-k-g]=rh;
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

	public static String UpdatePanGenomeAndMakeResultString(ClassificationResult cr, ConcurrentHashMap<Integer,Species> species)
	{
		String header = cr.read_header;
		String result = ",?,?|?";
		if (cr.species_kmer_hits!=null)
		{
			result = ","+cr.species_index+","+cr.species_kmer_hits.size()+"|"+cr.read_kmers;
			Species spp = species.get(cr.species_index);
			while(cr.species_kmer_hits.size()>0)
			{
				spp.ukmersfound.offer(cr.species_kmer_hits.pollFirst());
			}
			spp.readdepth++;
			species.put(cr.species_index,spp);
		}
		return header+result;
	}
	
	public static String[] parseSequenceFromFile(BufferedReader r) throws Exception
	{
		String[] out = null;
		String header = r.readLine(); 
		String fwd = r.readLine();
		if (fwd==null) return out;
		if (fwd.equals("")) return out;
		out = new String[2]; out[0]=header; out[1]=fwd;
		r.readLine();
		r.readLine();
		return out;
	}

	public static void main(String[] args) throws Exception
	{
		System.out.println();
		long timeZero = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		final int DEFAULT_BUFFER_SIZE=16384;

		String readfile="";
		String ukmdb = "";
		String outfilename="";
		int cores = Runtime.getRuntime().availableProcessors()-1;
		double pval_or_minfreq = 0.75d;
		int hll_log2m = 14;
		String line;
		
		for (int t=0; t<args.length; t++)
		{
			if (args[t].startsWith("s:")) pval_or_minfreq=Double.parseDouble(args[t].split(":")[1]);
			if (args[t].startsWith("o:")) outfilename=args[t].split(":")[1];
			if (args[t].startsWith("d:")) ukmdb=args[t].split(":")[1];
			if (args[t].startsWith("f:")) readfile=args[t].split(":")[1];
			if (args[t].startsWith("t:")) cores=Integer.parseInt(args[t].split(":")[1]);
			if (args[t].startsWith("l:")) hll_log2m=Integer.parseInt(args[t].split(":")[1]);
		}
		if (args==null || ukmdb.equals("") || readfile.equals("") || args[0].startsWith("-h") || args[0].startsWith("h"))
		{
			System.out.println("Please run the program as: \r\n \t java OCTOPUS d:database_folder f:fastq_file (can be gzipped) \r\n and optionally \r\n \t t:number_of_threads \r\n \t o:output_file_name \r\n \t s:probthreshold_or_minimum_hits (for classification, default is probability>0.75, any value >=1 will be minimum frequency of hits) \r\n \t l:log2m_value (for HyperLogLog) \r\n \t h or help or -h or -help to print this help");
			System.exit(0);
		}
		if (cores<1) {System.out.println("Cannot use less than one thread; will use one."); System.out.println(); cores = 1;}
		if (pval_or_minfreq<0) pval_or_minfreq=0.95d;
		if (hll_log2m<5 || hll_log2m>31) {System.out.println("Cannot use specified log_2 m value for HyperLogLog (either too large or too small); will try to use 14."); System.out.println(); hll_log2m = 14;}
		
		System.out.println("Max RAM allocated is "+allram/(1024*1024)+" MB"); System.out.println();
		System.out.println("Running program on "+cores+" threads");
		System.out.println();
		
		startTime = System.currentTimeMillis();
		System.out.println("Opening kmer:species database");
		KmerMetaDataBase store = new KmerMetaDataBase(ukmdb,pval_or_minfreq);
		int k = store.getKValue();
		boolean mask = store.getMaskValue();
		int alphabetSizePrime = 5;
		if (!mask) alphabetSizePrime = 4;
		final long [][] ppn = precomputedPowsForNucs(alphabetSizePrime,k);
		final HashMap<Character,Integer> nc = nucleotideCodes(mask);
		final HashMap<Character,Integer> ncr = nucleotideCodesRev(mask);
		long n = store.getNumKmers();
		System.out.print("\tThe database contains "+n+" k-mers ");
		Collection<Integer> keysc = store.getSpeciesKmerKeySet();
		ArrayList<Integer> keys = new ArrayList<Integer>(keysc);
		int nspecies = keys.size();
		System.out.println("and "+nspecies+" species");
		
		HyperLogLogPlus estimator = new HyperLogLogPlus(hll_log2m);
		long hllramused = (long)nspecies*estimator.getBytes().length;
		System.gc();
		usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		if ((double)(usedram+hllramused)/allram>0.8)
		{
			while((double)(usedram+hllramused)/allram>0.8)
			{
				hll_log2m = hll_log2m-1;
				estimator = new HyperLogLogPlus(hll_log2m);
				hllramused = (long)nspecies*estimator.getBytes().length;
			}
			System.out.println("\tProjected RAM usage is above max on-heap limit, setting log_2 m value for HyperLogLog to "+hll_log2m);
		}
		estimator=null;
		System.out.println("\tHyperLogLog log_2 m value set to "+hll_log2m);
		System.out.println("\t\tEstimated RAM usage of HyperLogLog counters for all species: "+hllramused/(1024*1024)+" MB.");
		System.out.println("\tClassification probability/frequency threshold value set to "+pval_or_minfreq);
		System.out.println();
		
		ConcurrentHashMap<Integer,Species> species = new ConcurrentHashMap(nspecies);
		for (int key : keys)
		{
			Species spp = new Species();
			spp.nkmers=store.getNumKmersSpecies(key);
			spp.ukmersfound=new HyperLogLogPlus(hll_log2m);
			spp.readdepth=0;
			species.put(key,spp);
		}
		
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("Database opened and species counter set up in "+elapsedTime/1000+" seconds");
		System.out.println();
		
		System.out.println("Parsing FASTQ read file and mapping genes");
		startTime = System.currentTimeMillis();
		String readOutFile = readfile.substring(0,readfile.lastIndexOf("."))+"_mappedReads.csv";
		if (!outfilename.equals("")) readOutFile = outfilename+"_mappedReads.csv";
		FileWriter rfilewriter = new FileWriter(readOutFile);
		BufferedWriter rwriter = new BufferedWriter(rfilewriter);
		rwriter.write("readId,taxonId,minimizKmersHitsTaxon|totMinimizKmers");
		rwriter.newLine();
		BufferedReader r = new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		if (readfile.endsWith(".gz")) {r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
		long i=0;
		long mapped=0;
		ConcurrentLinkedQueue<ClassificationResult> results = new ConcurrentLinkedQueue();
		boolean finished = false;
		Thread[] threadArray = new Thread[cores];
		for (int q=0; q<cores; q++)
		{
			String [] headerAndSequence = parseSequenceFromFile(r);
			if (headerAndSequence==null) {finished=true; break;}
			threadArray[q] = new Tentacle(q, results, headerAndSequence[0], headerAndSequence[1], k, mask, ppn, nc, ncr, store);
			threadArray[q].start();
			i++;
		}
		while (results.size()>0)
		{
			ClassificationResult cr = results.poll();
			String resultItem = UpdatePanGenomeAndMakeResultString(cr, species);
			if (resultItem.indexOf(",?")==-1) {mapped++;}
			rwriter.write(resultItem);rwriter.newLine();
		}
		boolean newThreadsStarted;
		while(!finished)
		{
			newThreadsStarted=false;
			for (int q=0; q<cores; q++)
			{
				if (!threadArray[q].isAlive())
				{
					String [] headerAndSequence = parseSequenceFromFile(r);
					if (headerAndSequence==null) {finished=true; break;}
					threadArray[q] = new Tentacle(q, results, headerAndSequence[0], headerAndSequence[1], k, mask, ppn, nc, ncr, store);
					threadArray[q].start();
					i++;
					newThreadsStarted=true;
				}
			}
			while (results.size()>0)
			{
				ClassificationResult cr = results.poll();
				String resultItem = UpdatePanGenomeAndMakeResultString(cr, species);
				if (resultItem.indexOf(",?")==-1) {mapped++;}
				rwriter.write(resultItem);rwriter.newLine();
			}
			if (i%6666==0 && newThreadsStarted)
			{
				System.gc();
				elapsedTime = System.currentTimeMillis() - startTime;
				System.out.println("\t"+i+" reads processed ("+mapped+" mapped) in "+elapsedTime/1000+" seconds");
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
				System.out.println("\t"+usedram/(1024*1024)+" MB RAM used ("+100*usedram/allram+"%)");
				System.out.println();
				for (int q=0; q<cores; q++) {threadArray[q].join();}
			}
		}
		for (int q=0; q<cores; q++) {threadArray[q].join();}
		while (results.size()>0)
		{
			ClassificationResult cr = results.poll();
			String resultItem = UpdatePanGenomeAndMakeResultString(cr, species);
			if (resultItem.indexOf(",?")==-1) {mapped++;}
			rwriter.write(resultItem);rwriter.newLine();
		}
		rwriter.flush();
		rwriter.close();
		r.close();
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("\t"+i+" reads processed ("+mapped+" mapped) in "+elapsedTime/1000+" seconds");
				System.gc();
		usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		System.out.println("\t"+usedram/(1024*1024)+" MB RAM used ("+100*usedram/allram+"%)");
		System.out.println("\t"+mapped+" reads assigned to taxa");
		System.out.println();
			
		System.out.println("Writing pan-genome classification and coverage summary on to file.. ");
		startTime = System.currentTimeMillis();
		System.out.print("\t");
		readOutFile = readfile.substring(0,readfile.lastIndexOf("."))+"_mappedGenomes.csv";
		if (!outfilename.equals("")) readOutFile = outfilename+"_mappedGenomes.csv";
		rfilewriter = new FileWriter(readOutFile);
		rwriter = new BufferedWriter(rfilewriter);
		rwriter.write("taxon_id,taxon_name,genome_coverage,read_depth");
		rwriter.newLine();
		keysc = species.keySet();
		keys = new ArrayList<Integer>(keysc);
		for (int key : keys)
		{
			Species spp = species.get(key);
			if (spp.readdepth>0)
			{
				rwriter.write(key+",");
				rwriter.write(store.getSpeciesName(key)+",");
				long cardi = spp.ukmersfound.cardinality();
				rwriter.write(cardi+"|"+spp.nkmers+",");
				rwriter.write(spp.readdepth+"");
				rwriter.newLine();
			}
		}			
		rwriter.flush();
		rwriter.close();
		store.closeKmerMetaDataBase();
		System.out.println("done.");
		System.out.println();
	
		elapsedTime = System.currentTimeMillis() - timeZero;
		System.out.println("Total time employed: "+elapsedTime/1000+" seconds");
	}

}

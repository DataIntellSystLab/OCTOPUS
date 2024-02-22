/*
del *.class
javac -cp ".;octopus_jars/*" GenomesToKmers.java
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\BVBRC_genomes.fasta.gz
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\who_priority_bacteria.fasta.gz
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\bacterial_refseq_urinary.fasta.gz
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\viral.1.1.genomic.fna.gz
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\mitochondrion.1.1.genomic.fna.gz
java -Xmx32G -cp ".;octopus_jars/*" GenomesToKmers ..\megares_database_v3.00.fasta.gz
*/


import java.io.*;
import java.io.File;
import java.lang.*;
import java.math.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.nio.*;
import java.nio.file.*;
import java.nio.charset.*;
import java.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.io.IOException;
import java.util.concurrent.TimeUnit;
import org.mapdb.*;
import org.mapdb.DB.*;
import org.mapdb.DBMaker.*;

public class GenomesToKmers
{
	public static void deleteFileOrFolder(File f) throws Exception
	{
		if (f.exists())
		{
			if (f.isDirectory())
			{
				File [] entries = f.listFiles();
				if (entries != null) {for (File entry : entries) {deleteFileOrFolder(entry);}}
			}
			f.delete();
		}
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
	public static long [] hashSequence(String fwd, int k, long [][] pp, HashMap<Character,Integer> nc, HashMap<Character,Integer> ncr, boolean mask)
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
		long[] minimizer_hashes;
		if ((fwd_hashes.length-step)>0)
		{
			minimizer_hashes = new long[fwd_hashes.length-step];
			long min = Long.MAX_VALUE;
			for (int g=0; g<minimizer_hashes.length; g++)
			{
				minimizer_hashes[g] = getMinimizer(fwd_hashes,min,g,step);
				min = minimizer_hashes[g];
			}
		}
		else
		{
			minimizer_hashes = new long[1];
			boolean found = false;
			long min = Long.MAX_VALUE;
			for (int g=0; g<fwd_hashes.length; g++)
			{
				if (fwd_hashes[g]>0 && fwd_hashes[g]<min) {min=fwd_hashes[g]; found=true;}
			}
			if (found) {minimizer_hashes[0] = min;} else  {minimizer_hashes[0] = -1;}
		}
		//System.out.println("minimize hashes = "+Arrays.toString(minimizer_hashes));
		return minimizer_hashes;
	}
	
	public static void updateKmersSpecies(String sequence, int k, ArrayList<Long> kmersSpecies, long [][] pp, HashMap<Character,Integer> nc, HashMap<Character,Integer> ncr, boolean mask)
	{
		long [] hashes = hashSequence(sequence, k, pp, nc, ncr, mask);
		for (int g=0; g<hashes.length; g++)
		{
			if (hashes[g]>=0) kmersSpecies.add(hashes[g]);
		}
	}
	
	public static void sortAndPrintArray(ArrayList<Long> kmersSpeciesL, int spp, String fileAndFolder) throws Exception
	{
		OutputStream out = new FileOutputStream(new File(fileAndFolder));
		BufferedWriter w = new BufferedWriter(new OutputStreamWriter(out, "UTF-8"),16384);
		if (kmersSpeciesL.size()>0)
		{
			Long [] kmersSpecies = kmersSpeciesL.toArray(new Long[kmersSpeciesL.size()]);
			Arrays.parallelSort(kmersSpecies);
			long prec = kmersSpecies[0];
			w.write(prec+","+spp); w.newLine();
			for (int i=1; i<kmersSpecies.length; i++)
			{
				long next = kmersSpecies[i];
				if (next!=prec)
				{
					prec=next;
					w.write(next+","+spp);
					w.newLine();
				}
			}
		}
		else {w.newLine();}
		w.flush();
		w.close();
	}

	public static long [] mergeSortedFile(String f1, String f2, String folder) throws Exception
	{
		BufferedReader r1 = new BufferedReader(new FileReader(folder+f1),16384);
		BufferedReader r2 = new BufferedReader(new FileReader(folder+f2),16384);
		long counts = 0;
		TimeUnit.MILLISECONDS.sleep(1);
		long outname = System.currentTimeMillis();
		BufferedWriter w = Files.newBufferedWriter(Paths.get(folder+outname), StandardCharsets.UTF_8);
		String line1 = r1.readLine();
		String line2 = r2.readLine();
		while (line1!=null && line2!=null)
		{
			if ( line1.equals("") || line2.equals("") || line1.equals("\n") || line2.equals("\n") || line1.equals("\r\n") || line2.equals("\r\n") ) {break;}
			long l1 = Long.parseLong(line1.split(",")[0]);
			long l2 = Long.parseLong(line2.split(",")[0]);
			if (l1<l2) {w.write(line1);w.newLine();line1=r1.readLine(); counts++;}
			else
			{
				if (l1>l2) {w.write(line2);w.newLine();line2=r2.readLine(); counts++;}
				else
				{
					if (l1==l2) 
					{
						String sppi1 = line1.split(",")[1];
						String sppi2 = line2.split(",")[1];
						String speciesList = sppi1+";"+sppi2;
						String [] sortedSpeciesListArray = speciesList.split(";");
						Arrays.sort(sortedSpeciesListArray);
						String sortedSpeciesList=sortedSpeciesListArray[0];
						for (int y=1; y<sortedSpeciesListArray.length; y++) {sortedSpeciesList+=";"+sortedSpeciesListArray[y];}
						if (!sppi1.equals(sppi2)) {w.write(l1+","+sortedSpeciesList); w.newLine(); counts++;}
						else {w.write(l1+","+sppi1); w.newLine(); counts++;}
						line1=r1.readLine();line2=r2.readLine();
					}
				}
			}
		}
		if ( line1==null || line1.equals("") || line1.equals("\n") || line1.equals("\r\n") ) {while(line2!=null){w.write(line2);w.newLine();line2=r2.readLine();counts++;}}
		if ( line2==null || line2.equals("") || line2.equals("\n") || line2.equals("\r\n") ) {while(line1!=null){w.write(line1);w.newLine();line1=r1.readLine();counts++;}}
		r1.close();
		r2.close();
		w.flush();
		w.close();
		deleteFileOrFolder(new File (folder+f1));
		deleteFileOrFolder(new File (folder+f2));
		long [] o = new long[2];
		o[0]=outname;
		o[1]=counts;
		return o;
	}

	public static void mergeKmerSpectraSpecies(String tempdir, String outf, int k, boolean mask, HashMap<Integer,String> speciesNames, int similthr) throws Exception
	{
		long startTime=System.currentTimeMillis();
		long elapsedTime=System.currentTimeMillis()-startTime;
		File f = new File(tempdir);
		long counts=0;
		String [] fileNames = f.list();
		LinkedList<String> f_even = new LinkedList();
		LinkedList<String> f_odd = new LinkedList();
		for (int i=0; i<fileNames.length; i++) {if (i%2==0) {f_even.add(fileNames[i]);}  else {f_odd.add(fileNames[i]);}}
		System.out.print("\t");
		while (f_odd.size()>0 && f_even.size()>0)
		{
			String f1 = f_odd.removeFirst();
			String f2 = f_even.removeFirst();
			long [] merge = mergeSortedFile(f1, f2, tempdir);
			String fm = merge[0]+"";
			if (f_odd.size()>f_even.size()) {f_even.add(fm);} else {f_odd.add(fm);}
			int siz=(f_odd.size()+f_even.size());
			if (siz%100==0)
			{
				System.out.print(siz+".. ");
				if (siz%1000==0)
				{
					System.out.println();
					elapsedTime=System.currentTimeMillis()-startTime;
					System.out.println("\t\t time elapsed: "+elapsedTime/1000+" seconds");
					System.out.print("\t");
				}
			}
			if (siz==1) counts = merge[1];
		}
		
		elapsedTime=System.currentTimeMillis()-startTime;
		System.out.println("done.");
		System.out.println("\tTime elapsed for merge: "+elapsedTime/1000+" seconds");
		String outfile = "";
		if (f_odd.size()>0) {outfile=f_odd.removeFirst();} else {outfile=f_even.removeFirst();}
		System.out.println();
		/*
		if (max_GB_usage>0)
		{
			double estimated_GB_usage = counts*0.000000008d;
			if (estimated_GB_usage>(double)max_GB_usage)
			{
				long newcounts = 0;
				System.out.println("\tRemoving random kmers to fit desired database size..");
				double prob_skip = 1d-(double)max_GB_usage/estimated_GB_usage;
				BufferedReader rc = new BufferedReader(new FileReader(tempdir+outfile),16384);
				String newOutFile = "reduced.txt";
				BufferedWriter w = Files.newBufferedWriter(Paths.get(tempdir+newOutFile), StandardCharsets.UTF_8);
				long linesRead=0;
				String l = rc.readLine();
				while(l!=null)
				{
					linesRead++;
					if (Math.random()>prob_skip) {w.write(l); w.newLine(); newcounts++;}
					l=rc.readLine();
					if (counts>=6) {if (linesRead%(counts/6)==0) {System.out.println("\t\t"+100*linesRead/counts+"%");}}
				}	
				rc.close();
				w.close();
				File df = new File(tempdir+outfile);
				df.delete();
				outfile = newOutFile;
				System.out.println("\t\tNumber of kmers reduced from "+counts+" to "+newcounts);
				counts = newcounts;
			}
			System.out.println();
		}
		*/
		boolean cluster=true;
		if (similthr<=0) {cluster=false;}
		HashMap<Integer,Integer> speciesToCluster = new HashMap<Integer,Integer>(speciesNames.size());
		TreeMap<Integer,String> clusterToHeader = new TreeMap<Integer,String>();
		if (cluster)
		{
			System.out.println("Clustering species at "+similthr+"% similarity");
			startTime=System.currentTimeMillis();
			HashMap<String,Integer> speciesCounts = new HashMap(speciesNames.size(),1.0f);
			HashMap<String,Integer> speciesCouplesMatches = new HashMap(3333333);
			DB db = DBMaker.fileDB(tempdir+"speciesCouplesMatchesDisk.db").fileMmapEnable().fileMmapEnableIfSupported().fileMmapPreclearDisable().allocateStartSize(666*1024*1024).allocateIncrement(666*1024*1024).cleanerHackEnable().make(); //.cleanerHackEnable(); db.getStore().fileLoad();
			HTreeMap<String,Integer> speciesCouplesMatchesDisk = db.hashMap("speciesCouplesMatchesDisk",Serializer.STRING,Serializer.INTEGER).counterEnable().create();
			boolean swapToDisk = false;
			BufferedReader rc = new BufferedReader(new FileReader(tempdir+outfile),16384);
			long linesRead=0;
			String l = rc.readLine();
			System.out.println("\tCalculating pairwise distances");
			while(l!=null)
			{
				float allram = (float)(Runtime.getRuntime().maxMemory());
				float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
				if (usedram/allram>0.8f){swapToDisk = true; System.out.println("\t\tMax RAM reached, swapping to disk.."); break;}
				linesRead++;
				String spp = l.split(",")[1];
				if (spp.indexOf(";")==-1)
				{
					Integer val = speciesCounts.get(spp);
					if (val!=null) {speciesCounts.put(spp,val+1);}
					else {speciesCounts.put(spp,1);}
				}
				else
				{
					String [] sppSet = spp.split(";");
					for (int y1=0; y1<sppSet.length-1; y1++)
					{
						Integer v1 = speciesCounts.get(sppSet[y1]);
						if (v1!=null) {speciesCounts.put(sppSet[y1],v1+1);}
						else {speciesCounts.put(sppSet[y1],1);}
						for (int y2=y1+1; y2<sppSet.length; y2++)
						{
							String pair = sppSet[y1]+";"+sppSet[y2];
							Integer v = speciesCouplesMatches.get(pair);
							if (v==null) {speciesCouplesMatches.put(pair,1);}
							else {speciesCouplesMatches.put(pair,v+1);}
						}
					}
					Integer val = speciesCounts.get(sppSet[sppSet.length-1]);
					if (val!=null) {speciesCounts.put(sppSet[sppSet.length-1],val+1);}
					else {speciesCounts.put(sppSet[sppSet.length-1],1);}
				}
				l=rc.readLine();
				if (counts>=6) {if (linesRead%(counts/6)==0){System.out.println("\t\t"+100*linesRead/counts+"% ("+usedram/1048576+" MB RAM used)");}}
			}
			if (swapToDisk)
			{
				while(l!=null)
				{
					linesRead++;
					String spp = l.split(",")[1];
					if (spp.indexOf(";")==-1)
					{
						Integer val = speciesCounts.get(spp);
						if (val!=null) {speciesCounts.put(spp,val+1);}
						else {speciesCounts.put(spp,1);}
					}
					else
					{
						String [] sppSet = spp.split(";");
						for (int y1=0; y1<sppSet.length-1; y1++)
						{
							Integer v1 = speciesCounts.get(sppSet[y1]);
							if (v1!=null) {speciesCounts.put(sppSet[y1],v1+1);}
							else {speciesCounts.put(sppSet[y1],1);}
							for (int y2=y1+1; y2<sppSet.length; y2++)
							{
								String pair = sppSet[y1]+";"+sppSet[y2];
								Integer v = speciesCouplesMatches.get(pair);
								if (v!=null) {speciesCouplesMatches.put(pair,v+1);}
								else
								{
									v = speciesCouplesMatchesDisk.get(pair);
									if (v!=null) {speciesCouplesMatchesDisk.put(pair,v+1);}
									else {speciesCouplesMatchesDisk.put(pair,1);}
								}
							}
						}
						Integer val = speciesCounts.get(sppSet[sppSet.length-1]);
						if (val!=null) {speciesCounts.put(sppSet[sppSet.length-1],val+1);}
						else {speciesCounts.put(sppSet[sppSet.length-1],1);}
					}
					l=rc.readLine();
					if (linesRead%666666==0) {db.commit();}
					if (counts>=6) {if (linesRead%(counts/6)==0){db.commit(); float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()); System.out.println("\t\t"+100*linesRead/counts+"% ("+usedram/1048576+" MB RAM used)");}}
				}
			}
			rc.close();
			int combinedSizes = speciesCouplesMatches.size()+speciesCouplesMatchesDisk.size();
			System.out.println("\t"+combinedSizes+" pairs of species with common kmers");
			ArrayList<String> toRemove = new ArrayList();
			for (Map.Entry<String,Integer> entry : speciesCouplesMatches.entrySet())
			{
				String pair = entry.getKey();
				String spp1 = pair.split(";")[0];
				String spp2 = pair.split(";")[1];
				int num = entry.getValue();
				int den = Math.min(speciesCounts.get(spp1),speciesCounts.get(spp2));
				int sim = (int)(100f*(float)num/(float)den);
				if (sim<=similthr) {toRemove.add(pair);}
			}
			for (int h=0; h<toRemove.size(); h++) {speciesCouplesMatches.remove(toRemove.get(h));}
           		toRemove = new ArrayList();
			for (Map.Entry<String,Integer> entry : speciesCouplesMatchesDisk.entrySet())
			{
				String pair = entry.getKey();
				String spp1 = pair.split(";")[0];
				String spp2 = pair.split(";")[1];
				int num = entry.getValue();
				int den = Math.min(speciesCounts.get(spp1),speciesCounts.get(spp2));
				int sim = (int)(100f*(float)num/(float)den);
				if (sim<=similthr) {toRemove.add(pair);}
			}
			for (int h=0; h<toRemove.size(); h++) {speciesCouplesMatchesDisk.remove(toRemove.get(h));}
           		combinedSizes = speciesCouplesMatches.size()+speciesCouplesMatchesDisk.size();
			System.out.println("\t"+combinedSizes+" pairs of species have similarity greater than "+similthr+"%");
			int oldspeciessize = speciesNames.size();
			if (speciesCouplesMatches.size()>0)
			{
				System.out.println("\tPerforming Chinese restaurant clustering");
				TreeSet<Integer> sppSet = new TreeSet<Integer>();
				Set<String> keys = speciesCouplesMatches.keySet();
				for (String key: keys)
				{
					int spp1 = Integer.parseInt(key.split(";")[0]);
					int spp2 = Integer.parseInt(key.split(";")[1]);
					sppSet.add(spp1);
					sppSet.add(spp2);
				}
				keys = speciesCouplesMatchesDisk.keySet();
				for (String key: keys)
				{
					int spp1 = Integer.parseInt(key.split(";")[0]);
					int spp2 = Integer.parseInt(key.split(";")[1]);
					sppSet.add(spp1);
					sppSet.add(spp2);
				}
				ArrayList<ArrayList<Integer>> clustering = new ArrayList();
				int p=0;
				System.out.print("\t\t");
				for (int e : sppSet)
				{
					boolean efound=false;
					p++;
					if (keys.size()>10) {if (p%(keys.size()/10)==0) {System.out.print(Math.round(100*(double)p/sppSet.size())+"% ");}}
					for (int h=0; h<clustering.size(); h++)
					{
						ArrayList currcluster=clustering.get(h);
						if (currcluster.contains(e)) {efound=true; break;}
						else 
						{
							int c=0;
							for (int ki=0; ki<currcluster.size(); ki++)
							{
								String c1=currcluster.get(ki)+";"+e;
								String c2=e+";"+currcluster.get(ki);
								if (speciesCouplesMatches.get(c1)!=null || speciesCouplesMatches.get(c2)!=null) {c++;}
							}
							if (c>=(currcluster.size()/2)) {currcluster.add(e); efound=true; break;}
						}
					}
					if (!efound)
					{
						ArrayList newcluster = new ArrayList();
						newcluster.add(e);
						clustering.add(newcluster);
					}
				}
				System.out.println(clustering.size()+" clusters found");
				Set<Integer> skeys = speciesNames.keySet();
				for (int g=0; g<clustering.size(); g++)
				{
					ArrayList<Integer> cclus = clustering.get(g);
					for (int gg=0; gg<cclus.size(); gg++)
					{
						int sppc = cclus.get(gg);
						speciesToCluster.put(sppc,g);
						String header = clusterToHeader.get(g);
						if (header==null) {clusterToHeader.put(g,speciesNames.get(sppc));}
						else {header=header+";"+speciesNames.get(sppc); clusterToHeader.put(g,header);}
						skeys.remove(sppc);
					}
				}
				int clcount=clustering.size();
				for (int key: skeys)
				{
					speciesToCluster.put(key,clcount);
					clusterToHeader.put(clcount,speciesNames.get(key));
					clcount++;
				}
			}
			if (clusterToHeader.size()>0)
				{System.out.println("\tFinal number of (clustered) species is "+clusterToHeader.size()+" (from "+oldspeciessize+")");}
			else {System.out.println("\tNo clusters created, keeping the number of species at "+oldspeciessize);}
			db.commit();
			db.close();
			Files.deleteIfExists(Paths.get(tempdir+"speciesCouplesMatchesDisk.db"));
			elapsedTime=System.currentTimeMillis()-startTime;
			System.out.println("\tTime employed for clustering: "+elapsedTime/1000+" seconds");
			System.out.println();
		}
		System.out.println("GZIPping final files..");
		startTime=System.currentTimeMillis();
		BufferedReader r = new BufferedReader(new FileReader(tempdir+outfile),16384);
		OutputStream zip = new GZIPOutputStream(new FileOutputStream(new File(outf+".gz")));
		BufferedWriter w = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"),16384);
		if (clusterToHeader.size()>0)
		{
			w.write("species="+clusterToHeader.size()); w.newLine();
			for (Map.Entry<Integer,String> entry : clusterToHeader.entrySet())
			{
				w.write(entry.getKey()+","+entry.getValue()); w.newLine();
			}
		}
		else
		{
			w.write("species="+speciesNames.size()); w.newLine();
			for (Map.Entry<Integer,String> entry : speciesNames.entrySet())
			{
				w.write(entry.getKey()+","+entry.getValue()); w.newLine();
			}
		}
		w.write("k="+(k+4)); w.newLine();
		w.write("mask_nonACGTchars="+mask); w.newLine();
		w.write("number_kmers="+counts); w.newLine();
		System.out.println("\t"+counts+" kmers in total");
		System.out.print("\t");
		String l = r.readLine();
		long linesWritten=0;
		while(l!=null)
		{
			if (clusterToHeader.size()>0)
			{
				String kmer = l.split(",")[0];
				String [] olds = l.split(",")[1].split(";");
				TreeSet specSet = new TreeSet();
				for (int ggg=0; ggg<olds.length; ggg++)
				{
					int news = speciesToCluster.get(Integer.parseInt(olds[ggg]));
					specSet.add(news);
				}
				String spec = ""+specSet.pollFirst();
				while(specSet.size()>0) {spec=spec+";"+specSet.pollFirst();}
				w.write(kmer+","+spec);
			}
			else {w.write(l);}
			w.newLine();
			l=r.readLine(); 
			linesWritten++;
			if (counts>=15) {if (linesWritten%(counts/15)==0) {System.out.print(100*linesWritten/counts+"% ");}}
		}
		r.close();
		w.flush();
		w.close();
		elapsedTime=System.currentTimeMillis()-startTime;
		System.out.print("done ("+elapsedTime/1000+" seconds). Deleting temp files..");
		File dir = new File(tempdir);
		for(File file: dir.listFiles()) if (!file.isDirectory()) file.delete();
		dir.delete();
		deleteFileOrFolder(new File(tempdir));
		System.out.println(" done.");
	}

	public static void main(String[] args) throws Exception
	{
		long timeZero = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;
		
		boolean mask = false;
		int k = 29;
		String dbfile="";
		int limit = 97;
		int similthr = 90;
		int samplefirstm = -1;
		
		if (args.length<=0) {System.out.println("\r\nType java GenomesToKmers -h for command options"); System.exit(0);}
		if (args.length==1)
		{	
			if (args[0].startsWith("-h"))
			{
				System.out.println("Please run the program as"); 
				System.out.println("\t(Windows) java -Xmx[desired_ram] -cp \".;octopus_jars/*\" GenomesToKmers genomes_fasta_file_or_folder (will use default parameters)");
				System.out.println("\t(UNIX/Linux) java -Xmx[desired_ram] -cp \".:octopus_jars/*\" GenomesToKmers genomes_fasta_file_or_folder (will use default parameters)");
				System.out.println("\t The input file can be also a fasta.gz");
				System.out.println("\t \t or alternatively run with options");
				System.out.println("\t f:genomes_fasta_file_or_folder ('f:' must be specified in this case) \r\n \t \t with one or more options among \r\n \t n:[y,n] (n:y if you want to mask non-ACGT characters with N, n:n otherwise that is the default) \r\n \t k:kmer_length (default 29, min 15, max 35 for ACGT and max 29 for ACGTN bases) \r\n \t s:percent_similarity (for clustering, default 90, <=0 for none) \r\n \t m:subsample_x (it selects only the first x species) \r\n \t -h or -help to print this help");
				System.exit(0);
			}
			else {dbfile=args[0];}
		}
		else
		{
			for (int t=0; t<args.length; t++)
			{
				if (args[t].startsWith("m:")) samplefirstm=Integer.parseInt(args[t].split(":")[1]);
				if (args[t].startsWith("f:")) dbfile=args[t].split(":")[1];
				if (args[t].startsWith("k:")) k=Integer.parseInt(args[t].split(":")[1]);
				if (args[t].startsWith("n:y")) mask=true;;
				if (args[t].startsWith("s:")) similthr=Integer.parseInt(args[t].split(":")[1]);
			}
		}
		if (dbfile.equals("")) {System.out.println("Please specify a genomes' input file in fasta format or a folder with all genomes \r\n \t (type java GenomesToKmers -h for command options)"); System.exit(0);}
		if (k%2==0) k=k+1;
		if (k<15) {System.out.println("Minimum value of k is 15"); k=15;} 
		if (mask && k>29) {System.out.println("Maximum value of k is 29"); k=29;}
		if (!mask && k>35) {System.out.println("Maximum value of k is 35"); k=35;}
		k = k-4;
		System.out.println();
		System.out.println("Length of kmers is set to "+(k+4)+" (minimizers are k-4, i.e. "+k+")");
		System.out.println();
		if (similthr>99) {System.out.println("Maximum value of similarity is 99%\r\n"); similthr=99;}
		if (!mask) System.out.println("Masking set to false: all kmers with non-ACGT character will be ignored");
		if (mask) System.out.println("Masking set to true: all kmers with non-ACGT character will use the additional base N");
		
		float allram = (float)(Runtime.getRuntime().maxMemory());
		System.out.println();
		System.out.println("Max RAM allocated is "+allram/1048576+" MB");
		System.out.println();
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		
		int alphabetSizePrime = 5;
		if (!mask) alphabetSizePrime = 4;
		long [][] pp = precomputedPowsForNucs(alphabetSizePrime,k);
		HashMap<Character,Integer> nc = nucleotideCodes(mask);
		HashMap<Character,Integer> ncr = nucleotideCodesRev(mask);
		
		String tempdir = "temp"+System.currentTimeMillis()+"/";
		new File(tempdir).mkdirs();
		
		String datafilerootname = dbfile.substring(Math.max(0,dbfile.lastIndexOf("/")+1));
		datafilerootname = dbfile.substring(Math.max(0,dbfile.lastIndexOf("\\")+1));
		int indperiod = datafilerootname.lastIndexOf(".");
		if (datafilerootname.indexOf(".gz")>-1 && datafilerootname.indexOf(".gz")<indperiod) indperiod=datafilerootname.indexOf(".gz");
		if (datafilerootname.indexOf(".fna")>-1 && datafilerootname.indexOf(".fna")<indperiod) indperiod=datafilerootname.indexOf(".fna");
		if (datafilerootname.indexOf(".fas")>-1 && datafilerootname.indexOf(".fas")<indperiod) indperiod=datafilerootname.indexOf(".fas");
		if (indperiod!=-1) {datafilerootname = datafilerootname.substring(0,indperiod);}	
		String ukf = datafilerootname+".KmerSpecies.txt";
		
		System.out.println("Reading species' database, creating kmer spectrum + mapping");
		startTime = System.currentTimeMillis();
		HashMap<Integer,String> speciesNames = new HashMap<Integer,String>();
		
		File dbfilefile = new File(dbfile);
		if (!dbfilefile.isDirectory())
		{
			BufferedReader r = new BufferedReader(new FileReader(dbfile),DEFAULT_BUFFER_SIZE);
			if(dbfile.endsWith(".gz")) {r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dbfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
			String header = r.readLine();
			while(!header.startsWith(">")) {header=r.readLine();}
			int i=0;
			long tempTime;
			long lapTime=0;
			while(true)
			{
				if (samplefirstm>0 && i>samplefirstm) {break;}
				tempTime=System.currentTimeMillis();
				if (header==null) break;
				if (!header.startsWith(">")) {System.out.println("Wrong FASTA format, exiting program."); System.exit(0);}
				String sequence=""; String lastkmer="";
				ArrayList<Long> kmersSpecies = new ArrayList();
				speciesNames.put(i,header.replaceAll(","," - "));
				int splits=0;
				while(true)
				{
					sequence = r.readLine();
					if (sequence==null) {header=null; break;}
					if (sequence.startsWith(">")) {header=sequence; break;}
					sequence=lastkmer+sequence;
					if (sequence.length()>=k) {updateKmersSpecies(sequence,k,kmersSpecies,pp,nc,ncr,mask);}
					if (sequence.length()>=k) {lastkmer=sequence.substring(sequence.length()-k,sequence.length());} else {lastkmer=sequence;}
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					if ( (100*usedram/allram)>limit )
					{
						System.out.println("\tOver "+limit+"% of available RAM used, printing results so far and resuming\r\n");
						splits++;
						sortAndPrintArray(kmersSpecies,i,tempdir+i+"_"+splits);
						kmersSpecies=new ArrayList();
						System.gc();
					}
				}
				sortAndPrintArray(kmersSpecies,i,tempdir+i);
				String mfn = i+"";
				for (int nsf=1; nsf<=splits; nsf++)
				{
					long [] merge = mergeSortedFile(mfn,i+"_"+nsf,tempdir);
					mfn = merge[0]+"";
				}
				i=i+1;
				lapTime+=System.currentTimeMillis()-tempTime;
				if (i%666==0 || (lapTime/1000)>1000)
				{
					System.gc();
					elapsedTime = System.currentTimeMillis() - startTime;
					System.out.println("\t"+(i+1)+" genomes processed in "+elapsedTime/1000+" seconds");
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					System.out.println("\t"+usedram/1048576+" MB RAM used ("+100*usedram/allram+"%)");
					System.out.println();
					lapTime=0;
				}
			}
			r.close();
		}
		if (dbfilefile.isDirectory())
		{
			String [] pathnames = dbfilefile.list();
			for (int i=0; i<pathnames.length; i++)
			{
				String species = dbfile+"/"+pathnames[i];
				BufferedReader rd = new BufferedReader(new FileReader(species),DEFAULT_BUFFER_SIZE);
				if(species.endsWith(".gz")) {rd=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(species),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
				String l = rd.readLine();
				if (l!=null)
				{
					while(!l.startsWith(">")) {l=rd.readLine();}
					if (!l.startsWith(">")) {System.out.println("Wrong FASTA format for file "+pathnames[i]+". Skipping it.");}
					else
					{
						ArrayList<Long> kmersSpecies = new ArrayList();
						speciesNames.put(i,l.replaceAll(","," - "));
						int splits=0;
						String lastkmer = "";
						while (true)
						{
							l=rd.readLine();
							if (l==null) {break;}
							if (!l.startsWith(">"))
							{
								l=lastkmer+l;
								if (l.length()>=k) {updateKmersSpecies(l,k,kmersSpecies,pp,nc,ncr,mask);}
								if (l.length()>=k) {lastkmer=l.substring(l.length()-k,l.length());} else {lastkmer=l;}
							}
							usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
							if ( (100*usedram/allram)>limit )
							{
								System.out.println("\tOver "+limit+"% of available RAM used, printing results so far and resuming\r\n");
								splits++;
								sortAndPrintArray(kmersSpecies,i,tempdir+i+"_"+splits);
								kmersSpecies=new ArrayList();
								System.gc();
							}
						}
						sortAndPrintArray(kmersSpecies,i,tempdir+i);
						String mfn = i+"";
						for (int nsf=1; nsf<=splits; nsf++)
						{
							long [] merge = mergeSortedFile(mfn,i+"_"+nsf,tempdir);
							mfn = merge[0]+"";
						}
					}
				}
				rd.close();
				if (i%(1+pathnames.length/100)==0)
				{
					System.gc();
					elapsedTime = System.currentTimeMillis() - startTime;
					System.out.println("\t"+(i+1)+" genomes processed in "+elapsedTime/1000+" seconds");
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					System.out.println("\t"+usedram/1048576+" MB RAM used ("+100*usedram/allram+"%)");
					System.out.println();
				}		
			}
		}
		
		System.gc();
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("\tAll genomes read and their kmers written in "+elapsedTime/1000+" seconds");
		usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		System.out.println("\t"+usedram/1048576+" MB RAM used ("+100*usedram/allram+"%)");
		System.out.println();
		
		System.out.println("Merging kmer spectra of all species");
		startTime = System.currentTimeMillis();
		mergeKmerSpectraSpecies(tempdir, ukf, k, mask, speciesNames, similthr);
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println();
		System.out.print("Time employed to print kmer spectrum for all linked species = "+elapsedTime/1000+" s\r\n");
		
		elapsedTime = System.currentTimeMillis() - timeZero;
		System.out.print("Total time employed  = "+elapsedTime/1000+" s\r\n");
	}
}

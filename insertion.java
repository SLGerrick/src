
import java.io.*;
public class insertion{
	String reference;
	String reference2;
	String donor;
	String donor2; //used when generating a new donor genome
	
	//constants/constraints
	int INSERT0=6950;
	int INSERT1=60000;
	int READ_LEN = 30;
	
	int readFrom=0;
	int perfect_match=0;
	int matched=0;
	int unmatched=0;
	int[] edges;
	int success_similar=0;
	int fail_similar=0;
	int success_perfect, success_sw=0;
	int fail_perfect, fail_sw=0;
	int found=0; 			//how many reads were insertions
	int found_right=0;		//how many reads were correctly identified as insertions
	int found_wrong=0;		//how many reads were incorrectly identified as insertions
	
	//for Smith Waterman algorithm
	int ind;
	double mu,delta;
	int GENOME_LENGTH=500000;
	int REF_SIZE=GENOME_LENGTH/5;
	int success=0;			//how many reads were successfully mapped
	int fail=0;				//how many reads were mapped incorrectly
	
	int[] runs= new int[5];
	int[] sw_insertion = new int[5]; 	//gives the index of any insertion found in a run through sw
	
	public double smithWaterman(String read, String ref_string, int run)
	{
	//	System.out.println("\nreadFrom: "+readFrom);
		mu=.5;
		delta=1;
		char[] ref=(reference.substring(run*REF_SIZE, (run+1)*REF_SIZE)+READ_LEN).toCharArray();
		char[] don=read.toCharArray();
		
		int ref_size = REF_SIZE;
		int read_size = READ_LEN;
		
		//score matrix
		double H[][]=new double[ref_size+1][read_size+1];
		
		//Index matrices
		double temp[]=new double[4];
		int I_i[][] = new int[ref_size+1][read_size+1];
		int I_j[][] = new int[ref_size+1][read_size+1];
		
		for(int i=1;i<=ref_size;i++){
		    for(int j=1;j<=read_size;j++){
		      temp[0] = H[i-1][j-1]+similarity_score(ref[i-1],don[j-1]); 
		      temp[1] = H[i-1][j]-delta;                  
		      temp[2] = H[i][j-1]-delta;                 
		      temp[3] = 0.;
		      H[i][j] = find_array_max(temp,4);
		      switch(ind){
		      case 0:                                  // score in (i,j) stems from a match/mismatch
		   	I_i[i][j] = i-1;
			I_j[i][j] = j-1;
			break;
		      case 1:                                  // score in (i,j) stems from a deletion in sequence A
		    I_i[i][j] = i-1;
			I_j[i][j] = j;
			break;
		      case 2:                                  // score in (i,j) stems from a deletion in sequence B
		      	I_i[i][j] = i;
			I_j[i][j] = j-1;
			break;
		      case 3:                                  // (i,j) is the beginning of a subsequence
		      	I_i[i][j] = i;
			I_j[i][j] = j;	
			break;
		      }
		    }
		  }

		  // find the best score
		  double H_max = 0.;
		  int i_max=0,j_max=0;
		  for(int i=1;i<=ref_size;i++){
		    for(int j=1;j<=read_size;j++){
		      if(H[i][j]>H_max){
				H_max = H[i][j];
				i_max = i;
				j_max = j;
		      }
		    }
		  }
		
		// Backtracking from H_max
		  int current_i=i_max;
		  int current_j=j_max;
		  int next_i=I_i[current_i][current_j];
		  int next_j=I_j[current_i][current_j];
		  int tick=0;
		  char consensus_a[] = new char[ref_size+read_size+2];
		  char consensus_b[]=new char[ref_size+read_size+2];

		  int insertion_index=-1;
		  Boolean has_insertion=false;
		  while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0))
		  {

			    if(next_i==current_i)
			    {
			    	consensus_a[tick] = '-';                  // insertion in reference genome
		//	    	System.out.println("Insertion");
			    	insertion_index=current_i;
			    	sw_insertion[run]=current_i;
			    	has_insertion=true;
			    }
			    	
			    else                   
			    	consensus_a[tick] = ref[current_i-1];   // match/mismatch in A
	
			    if(next_j==current_j) 
			    {
			    	consensus_b[tick] = '-';                  // deletion in B
			//    	System.out.println("Deletion");
			    	has_insertion=true;
			    }
			    else                   
			    	consensus_b[tick] = don[current_j-1];   // match/mismatch in B
	
			    current_i = next_i;
			    current_j = next_j;
			    next_i = I_i[current_i][current_j];
			    next_j = I_j[current_i][current_j];
			    tick++;
		    }
		  
		  int sw_index = i_max+REF_SIZE*run;
		  //if this is the correct bin
		  if(readFrom<REF_SIZE*(run+1) && readFrom>REF_SIZE*run)
		  {
			  //if the check is accurate
			  if(Math.abs(sw_index-readFrom)<READ_LEN+5)
			  {
				  //if it was an insertion
				  if((Math.abs(INSERT0-readFrom)<READ_LEN+5 && INSERT0-readFrom>0)||(Math.abs(INSERT1-readFrom)<READ_LEN+5 && INSERT1-readFrom>0)&&has_insertion)
				  {
					  found_right++;
				    	//System.out.println("Insertion found at: "+insertion_index);
				  }
				  success++;
			  }
			  else
			  {
				  fail++;
			  }
				  
		  }
		  runs[run]=sw_index;
		  return H_max;
	/*	  for(int i =0; i<32;i++)
		  {
			  System.out.print(consensus_a[i]);
		  }
		  System.out.println("");
		  for(int i =0; i<32;i++)
		  {
			  System.out.print(consensus_b[i]);
		  }
		  */
	}
	
	double similarity_score(char a,char b){

		  double result;
		  if(a==b){
		      result=1.;
		    }
		  else{
		      result=-mu;
		    }
		  return result;
	}
	
	double find_array_max(double array[],int length){

		  double max = array[0];            // start with max = first element
		  ind = 0;

		  for(int i = 1; i<length; i++){
		      if(array[i] > max){
			max = array[i];
			ind = i; 
		      }
		  }
		  return max;                    // return highest value in array
		}
	
	

	public insertion(File r, File d) {
		//TODO: get reference genome from a file
		reference2=ReadWriteTextFile.getContents(r);
		reference=reference2.substring(0, GENOME_LENGTH);
		donor2=generateDonor();
		/*		System.out.println(donor.length());
		try {
			ReadWriteTextFile.setContents(d, donor2);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
		donor = ReadWriteTextFile.getContents(d);
		donor = donor2.substring(0, GENOME_LENGTH);
		//donor2=generateDonor();
		 
		
	}
	
	//gets a read from the donor file
	public String getRead(){
		int x = (int)(Math.random()*GENOME_LENGTH);
		readFrom=x;
		
		if((Math.abs(INSERT0-readFrom)<READ_LEN+5 && INSERT0-readFrom>0)|| (Math.abs(INSERT1-readFrom)<READ_LEN+5&& INSERT1-readFrom>0))
		{
			found++;
	//		System.out.println("SEE ONE"+readFrom);
		}
		
		return donor.substring(x, x+READ_LEN);
	}
	
	
	//gets a read from the donor file at the specified index
	public String getRead(int i)
	{
		readFrom=i;
		if((Math.abs(INSERT0-readFrom)<READ_LEN+5 && INSERT0-readFrom>0)|| (Math.abs(INSERT1-readFrom)<READ_LEN+5&& INSERT1-readFrom>0))
		{
			found++;
	//		System.out.println("SEE ONE"+readFrom);
		}
		return donor.substring(i, i+READ_LEN);
	}
	
	public String generateDonor()
	{
		char[] buildDonor= reference2.toCharArray();
	
		int countSNP=0;
		
		//iterate through, and change the reference at each spot with probability 1.2%
		for(int i=0; i<reference2.length(); i++)
		{
			int rand = (int)(Math.random()*1000);
			if(rand<12)
			{
				if(rand<3)
					buildDonor[i]='a';
				else if(rand<6)
					buildDonor[i]='c';
				else if(rand<9)
					buildDonor[i]='t';
				else buildDonor[i]='g';
				countSNP++;
			}
		}
		System.out.println("SNPs: "+countSNP);
		String build = new String(buildDonor);
		
		//manual insertions
		String s0 = build.substring(0,INSERT0);
		String s1 = build.substring(INSERT0, INSERT1);
		String s2 = build.substring(INSERT1);
		build = s0+"a"+s1+"c"+s2;
		
		return build;
	}
	
	
	public int hasInsertion(String cur_read)
	{
		return hasInsertion(cur_read, 0);
	}
	
	
	public int hasInsertion(String cur_read, int index_initial)
	{
		String read0=cur_read.substring(0,READ_LEN/3);
		String read1=cur_read.substring(READ_LEN/3,2*READ_LEN/3);
		String read2=cur_read.substring(2*READ_LEN/3,READ_LEN);
		
		//find perfect matches
		Boolean succeeded=false;
		int index_perfect = reference.indexOf(cur_read);
		if(reference.contains(cur_read))
		{
			int index=0;
			perfect_match++;
			if(Math.abs(index_perfect-readFrom)<READ_LEN)
			{
				while(index_perfect>-1)
				{
					index=index_perfect;
					//System.out.println("perfect check "+readFrom);
					if(Math.abs(index_perfect-readFrom)<READ_LEN+5)
					{
						//System.out.println("perfect check found "+index);
						succeeded=true;
					}
					index_perfect = reference.indexOf(cur_read, index_perfect+1);
				}
			}
			if(succeeded)
				success_perfect++;
			else
			{
				fail_perfect++;
			//	System.out.println("PERFECT MATCH FAIL: "+index_perfect+" "+readFrom);
			}
			return index;
		}
		
		//imperfect matches
		else{
			int index = index_initial;
			int count0=0;
			int count1=0;
			int count2=0;
			
			//First 10 characters
			index = reference.indexOf(read0, index_initial);
			while(index>-1 && index<reference.length()-READ_LEN)
			{
				int match=0;
				String ref = reference.substring(index,index+READ_LEN);
				for(int i=0;i<10;i++)
				{
					if(ref.charAt(i+10)==read1.charAt(i))
						match++;
					if(ref.charAt(i+20)==read2.charAt(i))
						match++;
				}
				if(match>18)
				{
					matched++;
					if(Math.abs(index-readFrom)<READ_LEN)
						success_similar++;
					else fail_similar++;
					return index;
				}
				
				index=reference.indexOf(read0,index+1);
				count0++;
			}
			
			//Middle 10 characters	
			index = reference.indexOf(read1, index_initial);
			while(index>-1 && index<reference.length())
			{
				int match=0;
				String ref = reference.substring(index,index+READ_LEN);
				for(int i=0;i<10;i++)
				{
					if(ref.charAt(i)==read0.charAt(i))
						match++;
					if(ref.charAt(i+20)==read2.charAt(i))
						match++;
				}
				if(match>18)
				{
					matched++;
					if(Math.abs(index-readFrom)<READ_LEN)
						success_similar++;
					else fail_similar++;
					return index;
				}
				index=reference.indexOf(read1,index+1);
				count1++;
			}
			
			//Last 10 Characters
			index = reference.indexOf(read2, index_initial);
			while(index>-1 && index<reference.length())
			{
				int match=0;
				String ref = reference.substring(index,index+READ_LEN);
				for(int i=0;i<10;i++)
				{
					if(ref.charAt(i)==read0.charAt(i))
						match++;
					if(ref.charAt(i+10)==read1.charAt(i))
						match++;
				}
				if(match>18)
				{
					matched++;
					if(Math.abs(index-readFrom)<READ_LEN)
						success_similar++;
					else fail_similar++;
					return index;
				}
				index=reference.indexOf(read2,index+1);
				count2++;
			}
			
			//unmatched reads
			unmatched++;
			
			double best_score=-1;
			int best_run=-1;
			
			for(int i=0; i<5;i++)
			{
				double d = smithWaterman(cur_read, reference,i);
				if(d>best_score)
				{
					best_score=d;
					best_run=i;
					
				}
			}
			if(best_score==-1)
				System.out.println("TEST1\n\n");
			index=runs[best_run];
			int insertion_loc=sw_insertion[best_run];
			//check to see if it thinks it found an insertion
			if(insertion_loc > 0 && Math.abs(insertion_loc-INSERT0)>4 && Math.abs(insertion_loc-INSERT0)>4)
			{
				System.out.println("wrong insertion found at: "+insertion_loc);
				found_wrong++;
			}
			
			if(Math.abs(index-readFrom)<READ_LEN)
				success_sw++;
			else fail_sw++;
			return index;
			
			//System.out.println("---"+readFrom);
		}
	}
	
	

	public static void main(String[] args){
		//TODO read a random set of characters from the donor genome
		System.out.println("Start");
		
		
		File ref_file= new File("/Users/slgerrick/Desktop/CS124/ref");//genome file
		File donor_file= new File("/Users/slgerrick/Desktop/CS124/refcopy2");

		insertion g=new insertion(ref_file, donor_file);
		int num_reads=100;
		
		long start = System.currentTimeMillis();
		
		
		for(int i=0; i<num_reads;i++)
		{
			String m_read = g.getRead();
			int a = g.hasInsertion(m_read);
			
			//Baseline
			//int x=g.INSERT0-g.READ_LEN+i-1;
		/*	String m_read = g.getRead(x);
			
			System.out.println("read: "+i);
			int best_run=-1;
			double best_score=-1;
			for(int j=0;j<5;j++)
			{
				//System.out.println("check:"+j);
				double d=g.smithWaterman(m_read, g.reference, j);
				if(d>best_score)
				{
					best_score=d;
					best_run=j;
					
				}
			}
			int a=g.runs[best_run];
			*/
			//checks for accuracy
			  if(Math.abs(a-g.readFrom)<g.READ_LEN+5)
			  {
				//check if an insertion was found
				  if((Math.abs(g.INSERT0-a)<(2*g.READ_LEN) && g.INSERT0-a>0) || (Math.abs(g.INSERT1-a)<(2*g.READ_LEN) && g.INSERT1-a>0))
				  {
					  g.found_right++;
				  }
			  }
			  

		}
		System.out.println("Found right: "+g.found_right);
		
		long end = System.currentTimeMillis();
		long time = end-start;
		
		long start2 = System.currentTimeMillis();
		
		System.out.println("Perfect: \t"+g.perfect_match+" Matched: "+g.matched+" Unmatched: "+g.unmatched);
		System.out.println("\tSW all\tSW in\tSim\t Perfect");
		System.out.println("Succeeded:\t"+g.success_sw+"\t"+g.success+"\t"+g.success_similar+"\t"+g.success_perfect+" \nFailed: \t"+g.fail_sw+"\t"+g.fail+"\t"+g.fail_similar+ "\t"+g.fail_perfect+" \nFound: \t"+g.found_right +"/ "+g.found+" with :"+g.found_wrong+" misidentified");
		System.out.println(""+time+" Milliseconds to perform"+num_reads+" reads");
		
		
		for(int i=0; i<g.READ_LEN+2; i++)
		{
			int x=g.INSERT0-g.READ_LEN+i-1;
			String m_read=g.getRead(x);
			int a=g.hasInsertion(m_read);
			System.out.println("at: "+a);
		}
	
		System.out.println("Perfect: \t"+g.perfect_match+" Matched: "+g.matched+" Unmatched: "+g.unmatched);
		System.out.println("\tSW all\tSW in\tSim\t Perfect");
		System.out.println("Succeeded:\t"+g.success_sw+"\t"+g.success+"\t"+g.success_similar+"\t"+g.success_perfect+" \nFailed: \t"+g.fail_sw+"\t"+g.fail+"\t"+g.fail_similar+ "\t"+g.fail_perfect+" \nFound: \t"+g.found_right +"/ "+g.found+" with :"+g.found_wrong+" misidentified");

		
		end = System.currentTimeMillis();
		time = end-start2;
		System.out.println(""+time+" Milliseconds");
		
		
/*		
		int success_check=0;
		for(int i=0; i<READ_LEN; i++)
		{
			int x=g.INSERT1-READ_LEN+i+1;
			String m_read=g.getRead(x);
			int a=g.hasInsertion(m_read);
			
			if(Math.abs(g.INSERT0-a)<READ_LEN+5)
			{
				success_check++;
				//a=g.hasInsertion(m_read, a+1);
			}
		}
		System.out.println("Success Check: "+success_check+" ");
		
		System.out.println("Perfect: \t"+g.perfect_match+" Matched: "+g.matched+" Unmatched: "+g.unmatched);
		
		System.out.println("\t\tSW all\tSW in\tSim\t Perfect");
		System.out.println("Succeeded:\t"+g.success_sw+"\t"+g.success+"\t"+g.success_similar+"\t"+g.success_perfect+" \nFailed: \t"+g.fail_sw+"\t"+g.fail+"\t"+g.fail_similar+ "\t"+g.fail_perfect+" \nFound: \t"+g.found);
		*/
		
	/*	int success_check=0;
		for(int i=0; i<READ_LEN; i++)
		{
			int x=g.INSERT1-READ_LEN+i+1;
			String m_read=g.getRead(x);
			int a=g.hasInsertion(m_read);
			while(a!=-1)
			{
				if(Math.abs(g.INSERT0-a)<READ_LEN+5)
				{
					success_check++;
					a=g.hasInsertion(m_read, a+1);
				}
			}
		}
		System.out.println("Success Check: "+success_check);
		*/
		
		
		System.out.println("\nDone");
		
	}
	
		

}

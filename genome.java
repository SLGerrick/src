
public class genome {
	String reference;
	int[] edges;
	static int NUMBINS=6667;
	
	public genome(String f){
		//TODO: get reference genome from a file
		reference=f;
	}
	
	
	//takes in a read of length 30 and returns the most likely bin that this read is stored in.  If no bins match, it returns -1
	public int getBin(String read30){
		
		//TODO read and map to the file, if the match is perfect
		String read0=read30.substring(0, 10);
		String read1=read30.substring(10, 20);
		String read2=read30.substring(20, 30);
		String cur_ref;
		
		//for 15 portions
		int places1[]=new int[50];
		int count1 =0;
		int places2[]=new int[50];
		int count2 =0;
		
		for(int i=0;i<reference.length()-30;i++){
			cur_ref=reference.substring(i, i+30);
			if(read30.substring(0,30).equals(cur_ref))
				return i/NUMBINS;
			if(cur_ref.contains(read0) & cur_ref.contains(read1))
				return i/NUMBINS;
			if(cur_ref.contains(read0) & cur_ref.contains(read2))
				return i/NUMBINS;
			if(cur_ref.contains(read1) & cur_ref.contains(read2))
				return i/NUMBINS;
			
		}
		
		//if read30 is not matched, check for a match on either half (running concurrently), if half matches, then find where on the other half the match stops working
		
		//assumes that a read of length 15 will not match to more than 15 locations

		for(int i=0;i<reference.length()-30;i++){
			//if half the read is a match, see how far into the second half the match goes
			if(read30.substring(0,15).equals(reference.substring(i, i+15)))
			{
				for(int k=0; k<15;k++)
				{
					if(read30.charAt(15+k)!=reference.charAt(i+k))
					{
						places1[count1]=i+15+k;
						count1++;
						break;
					}
				}
					
			}
			if(read30.substring(0,15).equals(reference.substring(i, i+15)))
			{
				for(int k=0; k<15;k++)
				{
					if(read30.charAt(15-k)!=reference.charAt(i+k))
					{
						places2[count2]=i+15-k;
						count2++;
						break;
					}
				}
			}
		}
		
		//determine which place the change occurs in
		
		return -1;
	}
	/*
	public static void main(String[] args){
		//TODO read a random set of characters from the donor genome
		String m_read="actactatctagatagacatagca";
		String f = "";//genome file
		
		genome g=new genome(f);
		
		int[] bin = new int[NUMBINS];
		int i = g.getBin(m_read);
		if(i>-1)
			bin[i]++;
		else
		{
			int j=g.getBin(m_read.substring(0, 15));
			int k=g.getBin(m_read.substring(15));
			if(j>-1&& k==-1)
				bin[j]++;
			else if(j==-1 && k>-1)
				bin[k]++;
			else
				System.out.println("ERROR\n");
		}
	}
	
		*/

}

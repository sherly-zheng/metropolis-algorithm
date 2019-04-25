import java.util.concurrent.ThreadLocalRandom;

public class MetropolisThread extends Thread{
	private int N_m;
	private int N_f;
	private int n;
	private double B;
	private double C;
	private double T;
	private int rand_index;
	private double energy_diff;
	private int[] config;
	private double[] m;
	private double[] cp;
	private double m_mean;
	private double cp_mean;
	
	public MetropolisThread(int num_metro, int num_flips, int num_spins, double b_param, double c_param, double temperature){
			N_m = num_metro;
			N_f = num_flips;
			n = num_spins;
			B = b_param;
			C = c_param;
			T = temperature;
			rand_index = 0;
			energy_diff = 0;
			config = new int[n];	//	current configuration
			m = new double[N_m];	// 	stores the magnetization of each optimal configuration. We have N_m optimal configurations since we make N_m calls to MetropolisAlgorithm().
			cp = new double[N_m];	// 	stores the pair correlation of each optimal configuration. We have N_m optimal configurations since we make N_m calls to MetropolisAlgorithm().
			m_mean = 0;				//	the mean magnetization for this thread.
			cp_mean = 0;			//	the mean pair correlation for this thread.
	}
	
	// calls the MetropolisAlgorithm() N_m times, and computes mean magnetization and mean pair correlation for this thread.
	public void run(){		
		double m_sum = 0;
		double cp_sum = 0;
		for(int i = 0; i < N_m; i++){
			MetropolisAlgorithm();
			computeMagAndPairCorr(i);
			m_sum += m[i];
			cp_sum += cp[i];
		}
		m_mean = m_sum/N_m;
		cp_mean = cp_sum/N_m;
	}
		
	public double getMagMean(){
		return m_mean;
	}
	
	public double getPairCorrMean(){
		return cp_mean;
	}
	
	/*	computes magnetization and pair correlation of the optimal configuration which is the result of the i^th call to MetropolisAlgorithm(). 
		Stores the resulting magnetization in m[i] and pair correlation in cp[i].
	*/
	public void computeMagAndPairCorr(int i){
		int m_sum = 0;
		int cp_sum = 0;
		
		for(int j = 0; j < n; j++){
			m_sum += config[j];
			cp_sum += config[j]*config[getNextIndexOf(j)];
		}
		m[i] = (m_sum*1.0)/n;
		cp[i] = (cp_sum*1.0)/n;
	}

	//	performs the Metropolis Algorithm. The result of this function is the "most optimal" configuration.
	public void MetropolisAlgorithm(){
		initConfig();
		boolean replaceConfig;
		int flip_count = 0;
		while(flip_count < (n * N_f)){
			selectRandomSpin();
			computeEnergyDiff();
			replaceConfig = canReplaceConfig();
			if(replaceConfig)
				config[rand_index] = -config[rand_index];
			flip_count++;
		}
	}
	
	public void initConfig(){
		int value;
		if(C >= 0)
			value = 1;
		else 
			value = -1;
		for(int i = 0; i < n; i++)
			config[i] = (int)Math.pow(value, i);
	}

	public void selectRandomSpin(){
		rand_index = ThreadLocalRandom.current().nextInt(0, n);
	}
	
	/*	returns the index that comes before index i. If i == 0, its "previous" index will be the last index, aka n-1, in the config array. 
		In the case where i == 0, (i-1)%n will give us -1 instead of n-1 because Java treats % as a remainder operator instead of a
		modulus operator. By returning ((i-1) + n) % n, we are guaranteed a positive integer index for all i.
	*/
	public int getPrevIndexOf(int i){
		return ((i-1) + n) % n;
	}
	
	//	returns the index that comes after index i. If i == the last index in the config array, its "next" index will be 0.
	public int getNextIndexOf(int i){
		return (i+1) % n;
	}
	
	/*	computes the energy difference between current configuration and flipped configuration.
		Let's call the current configuration sigma_0.
		In sigma_0, we select one spin at random to be flipped. Let's call that spin s_i. 
		The configuration resulting from that flip will be called sigma_1.
		
		The energy of a configuration, E(sigma) = -((B*s_1 + C*s_1*s_2) + (B*s_2 + C*s_2*s_3) + ... + (B*s_n + C*s_n*s_n+1)).
		So, 
			E(sigma_1)	= 	-((B*s_1 + C*s_1*s_2) + ... + (B*s_i-1 + C*s_i-1*(-s_i)) + (B*(-s_i) + C*(-s_i)*s_i+1) + ... + (B*s_n + C*s_n*s_n+1)) 
		and	
			E(sigma_0)	= 	-((B*s_1 + C*s_1*s_2) + ... + (B*s_i-1 + C*s_i-1*  s_i ) + (B*  s_i  + C*  s_i *s_i+1) + ... + (B*s_n + C*s_n*s_n+1)).
		
		We obtain the energy difference between the two configurations by calculating E(sigma_1) - E(sigma_0).
			E(sigma_1) - E(sigma_0) =  	-((B*s_1 + C*s_1*s_2) + ... + (B*s_i-1 + C*s_i-1*(-s_i)) + (B*(-s_i) + C*(-s_i)*s_i+1) + ... + (B*s_n + C*s_n*s_n+1)) 
										-
										-((B*s_1 + C*s_1*s_2) + ... + (B*s_i-1 + C*s_i-1*  s_i ) + (B*  s_i  + C*  s_i *s_i+1) + ... + (B*s_n + C*s_n*s_n+1)).
										
		Since none of the spins except s_i change/flip, subtracting E(sigma_0) from E(sigma_1) will result in all summands being canceled out, EXCEPT for those where we see s_i. 
		The result is,
			E(sigma_1) - E(sigma_0) =	- ((B*s_i-1 + C*s_i-1*(-s_i)) + (B*(-s_i) + C*(-s_i)*s_i+1)) 
										+ ((B*s_i-1 + C*s_i-1*  s_i ) + (B*  s_i  + C*  s_i *s_i+1)).
		
		Flip it around to get,
			E(sigma_1) - E(sigma_0) =	  ((B*s_i-1 + C*s_i-1*  s_i ) + (B*  s_i  + C*  s_i *s_i+1)) 
										- ((B*s_i-1 + C*s_i-1*(-s_i)) + (B*(-s_i) + C*(-s_i)*s_i+1))
		
									=	   (B*s_i-1 + C*s_i-1*  s_i ) + (B*  s_i  + C*  s_i *s_i+1)
										-  (B*s_i-1 - C*s_i-1*  s_i ) - (-B* s_i  - C*  s_i *s_i+1)).
		
		The values for w, x, y and z in the code below are equivalent to,
			w = B*s_i-1
			x = C*s_i-1*s_i
			y = B*s_i
			z = C*s_i *s_i+1	
	*/
	public void computeEnergyDiff(){
		int s_i = config[rand_index];
		int s_i_prev = config[getPrevIndexOf(rand_index)];
		int s_i_next = config[getNextIndexOf(rand_index)];
		
		double w = B * s_i_prev;
		double x = C * s_i_prev * s_i;
		double y = B * s_i;
		double z = C * s_i * s_i_next;
		energy_diff = (w+x) + (y+z) - (w-x) - (-y-z);
	}
	
	//	tests if the flipped configuration can replace current configuration.
	public boolean canReplaceConfig(){
		if(energy_diff <= 0)
			return true;
		else{
			double p = Math.exp((-energy_diff) / T);
			double r = ThreadLocalRandom.current().nextDouble(0, 1.0);
			if(r < p)
				return true;
			return false;
		}
	}
}
	
import java.util.Scanner;
import java.io.*;
public class Main{
	public static void main(String args[]){
		try{
			//	get user inputs
			int n;
			double B, C, T;
			String file_name;	
			Scanner scanner = new Scanner(System.in);
			
			System.out.print("Enter the output text file name: ");
			file_name = scanner.next();	
			BufferedWriter outfile = new BufferedWriter(new FileWriter(file_name));
			System.out.print("Enter a value for n: ");
			n = scanner.nextInt();
			System.out.print("Enter a value for B: ");
			B = scanner.nextDouble();
			System.out.print("Enter a value for C: ");
			C = scanner.nextDouble();
			System.out.print("Enter a value for T: ");
			T = scanner.nextDouble();
				
			//	threading
			int N_t = 1000;
			int N_m = 90;
			int N_f = 35;
			MetropolisThread[] thread = new MetropolisThread[N_t];
			
			for(int i = 0; i < N_t; i++){
				thread[i] = new MetropolisThread(N_m, N_f, n, B, C, T);
				thread[i].start();
			}
			
			for(int i = 0; i < N_t; i++){
				try{
					if(thread[i].isAlive())
						thread[i].join();
				}
				catch (Exception e){
					System.out.println(e.getMessage());
				}
			}
			
			//	compute global magnetization mean and global pair correlation mean.
			double mu_m_sum = 0;
			double mu_cp_sum = 0;
			double mu_m = 0;						//	global magnetization mean.
			double mu_c = 0;						//	global pair correlation mean.
			double[] m = new double[N_t];			//	stores magnetization means of each thread.
			double[] c = new double [N_t];			//	stores pair correlation means of each thread.
 			
			for(int j = 0; j < N_t; j++){
				m[j] = thread[j].getMagMean();
				c[j] = thread[j].getPairCorrMean();
				mu_m_sum += m[j];
				mu_cp_sum += c[j];
				
				//	output m_means[] and cp_means[] to file for statistical error analysis.
				outfile.write("m[" + j + "] = " + m[j] + "     c[" + j + "] = " + c[j]);
				outfile.newLine();
			} 
			mu_m = mu_m_sum/N_t;
			mu_c = mu_cp_sum/N_t;
			outfile.write("mu_m = " + mu_m + "     mu_c = " + mu_c);
			
			scanner.close();
			outfile.close();
		}
		catch(IOException e){
			System.out.println(e.getMessage());
		}
	}
}
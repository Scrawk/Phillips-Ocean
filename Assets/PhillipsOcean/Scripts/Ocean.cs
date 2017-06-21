using UnityEngine;
using System;
using System.Collections;
using System.Threading;

namespace PhillipsOcean
{
	/// <summary>
	/// The ocean game object.
	/// Is responsible for creating the ocean mesh and other
	/// resources.
	/// </summary>
	public class Ocean : MonoBehaviour 
	{	

		const float GRAVITY = 9.81f;

		/// <summary>
		/// The material used to render the mesh. The actual rendering is very basic. 
		/// This project is more about how to generate the spectrum.
		/// </summary>
		[SerializeField]
		Material m_mat;

		/// <summary>
		/// The number of meshes on the x and z axis.
		/// </summary>
		[SerializeField]
		int m_numGridsX = 2;
		[SerializeField]
		int m_numGridsZ = 2;

		/// <summary>
		/// The fourier size. Must be a pow2 number.
		/// </summary>
		[SerializeField]
		int N = 64;

		/// <summary>
		/// The world space length of each mesh.
		/// </summary>
		[SerializeField]
		float m_length = 64;

		/// <summary>
		/// The phillips spectrum parameter. Affects heights of waves.
		/// Setting this to hight can cause the waves to curl back
		/// on themselves.
		/// </summary>
		[SerializeField]
		float m_waveAmp = 0.0002f;

		/// <summary>
		/// The wind speed and the direction.
		/// </summary>
		[SerializeField]
		Vector2 m_windSpeed = new Vector2(32.0f,32.0f);

		/// <summary>
		/// The game object holding each mesh.
		/// </summary>
		GameObject[] m_oceanGrid;

		/// <summary>
		/// The actual mesh. There is only one mesh generated.
		/// The mesh is then tiled in the x and z axis.
		/// </summary>
		Mesh m_mesh;

		/// <summary>
		/// Just the fourier size plus 1.
		/// </summary>
		int Nplus1;	

		/// <summary>
		/// The wind direction.
		/// </summary>
		Vector2 m_windDirection;

		/// <summary>
		/// Use to perform the fourier transform.
		/// </summary>
		FourierCPU m_fourier;

		/// <summary>
		/// Buffers for fourier transform.
		/// </summary>
		Vector2[,] m_heightBuffer;
		Vector4[,] m_slopeBuffer, m_displacementBuffer;

		/// <summary>
		/// Holds the spectrum which will be then be transformed into
		/// the wave heights/slopes.
		/// </summary>
		Vector2[] m_spectrum, m_spectrum_conj;

		/// <summary>
		/// 
		/// </summary>
		Vector3[] m_position;
		Vector3[] m_vertices;
		Vector3[] m_normals;

		/// <summary>
		/// Just holds so data for the spectrum that can be precomputed.
		/// </summary>
		float[] m_dispersionTable;

		/// <summary>
		/// The fresnel look up table.
		/// </summary>
		Texture2D m_fresnelLookUp;

		/// <summary>
		/// Has the thread generating the data finished.
		/// </summary>
		volatile bool done = true;

		void Start() 
		{
			Nplus1 = N+1;
			
			m_fourier = new FourierCPU(N);
			
			m_windDirection = new Vector2(m_windSpeed.x, m_windSpeed.y);
			m_windDirection.Normalize();
			
			m_dispersionTable = new float[Nplus1*Nplus1];
			
			for (int m_prime = 0; m_prime < Nplus1; m_prime++) 
			{
				for (int n_prime = 0; n_prime < Nplus1; n_prime++) 
				{
					int index = m_prime * Nplus1 + n_prime;
					m_dispersionTable[index] = Dispersion(n_prime,m_prime);
				}
			}
			
			m_heightBuffer = new Vector2[2,N*N];
			m_slopeBuffer = new Vector4[2,N*N];
			m_displacementBuffer = new Vector4[2,N*N];
			
			m_spectrum = new Vector2[Nplus1*Nplus1];
			m_spectrum_conj = new Vector2[Nplus1*Nplus1];

			m_position = new Vector3[Nplus1*Nplus1];
			m_vertices = new Vector3[Nplus1*Nplus1];
			m_normals = new Vector3[Nplus1*Nplus1];
			
			m_mesh = MakeMesh(Nplus1);
			;
			
			m_oceanGrid = new GameObject[m_numGridsX*m_numGridsZ];

			for(int x = 0; x < m_numGridsX; x++)
			{
				for(int z = 0; z < m_numGridsZ; z++)
				{
					int idx = x + z * m_numGridsX;
					
					m_oceanGrid[idx] = new GameObject("Ocean grid " + idx.ToString());
					m_oceanGrid[idx].AddComponent<MeshFilter>();
					m_oceanGrid[idx].AddComponent<MeshRenderer>();
					m_oceanGrid[idx].GetComponent<Renderer>().material = m_mat;
					m_oceanGrid[idx].GetComponent<MeshFilter>().mesh = m_mesh;
					m_oceanGrid[idx].transform.Translate(new Vector3(x * m_length - m_numGridsX*m_length/2, 0.0f, z * m_length - m_numGridsZ*m_length/2));
					m_oceanGrid[idx].transform.parent = this.transform;

				}
			}

            UnityEngine.Random.InitState(0);
			
			Vector3[] vertices = m_mesh.vertices;
			
			for (int m_prime = 0; m_prime < Nplus1; m_prime++) 
			{
				for (int n_prime = 0; n_prime < Nplus1; n_prime++) 
				{
					int index = m_prime * Nplus1 + n_prime;
		
					m_spectrum[index] = GetSpectrum( n_prime,  m_prime);
					
					m_spectrum_conj[index] = GetSpectrum(-n_prime, -m_prime);
					m_spectrum_conj[index].y *= -1.0f;
		
					m_position[index].x = vertices[index].x =  n_prime * m_length/N;
					m_position[index].y = vertices[index].y =  0.0f;
					m_position[index].z = vertices[index].z =  m_prime * m_length/N;

				}
			}
			
			m_mesh.vertices = vertices;
			m_mesh.RecalculateBounds();
			
			CreateFresnelLookUp();
		}

		/// <summary>
		/// Create a fresnel lookup table. This is the formula
		/// to calculate a materials fresnel value based on 
		/// its refractive index. Since its a little math heavy
		/// a look up table is used rather than caculate it in 
		/// the shader. In practise this method does not look any better
		/// than cheaper approximations but is included out of interest.
		/// </summary>
		void CreateFresnelLookUp()
		{
			float nSnell = 1.34f; //Refractive index of water
		
			m_fresnelLookUp = new Texture2D(512, 1, TextureFormat.Alpha8, false);
			m_fresnelLookUp.filterMode = FilterMode.Bilinear;
			m_fresnelLookUp.wrapMode = TextureWrapMode.Clamp;
			m_fresnelLookUp.anisoLevel = 0;
			
			for(int x = 0; x < 512; x++)
			{
				float fresnel = 0.0f;
				float costhetai = (float)x/511.0f;
				float thetai = Mathf.Acos(costhetai);
				float sinthetat = Mathf.Sin(thetai)/nSnell;
				float thetat = Mathf.Asin(sinthetat);
				
				if(thetai == 0.0f)
				{
					fresnel = (nSnell - 1.0f)/(nSnell + 1.0f);
					fresnel = fresnel * fresnel;
				}
				else
				{
					float fs = Mathf.Sin(thetat - thetai) / Mathf.Sin(thetat + thetai);
					float ts = Mathf.Tan(thetat - thetai) / Mathf.Tan(thetat + thetai);
					fresnel = 0.5f * ( fs*fs + ts*ts );
				}
				
				m_fresnelLookUp.SetPixel(x, 0, new Color(fresnel,fresnel,fresnel,fresnel));
			}
			
			m_fresnelLookUp.Apply();
			
			m_mat.SetTexture("_FresnelLookUp", m_fresnelLookUp);
		}
		
		/// <summary>
		/// Evaulate the waves for time period each update.
		/// </summary>
		void Update () 
		{

			//If still running return.
			if(!done) return;

			//Set data generated form last calculations.
			m_mesh.vertices = m_vertices;
			m_mesh.normals = m_normals;
			m_mesh.RecalculateBounds();

			//Start new calculations for time period t.
			done = false;
			
			Nullable<float> time = Time.realtimeSinceStartup;
			
			ThreadPool.QueueUserWorkItem(new WaitCallback(RunThreaded), time);
		}

		/// <summary>
		/// Gets the spectrum vaule for grid position n,m. 
		/// </summary>
		Vector2 GetSpectrum(int n_prime, int m_prime) 
		{
			Vector2 r = GaussianRandomVariable();
			return r * Mathf.Sqrt(PhillipsSpectrum(n_prime, m_prime) / 2.0f);
		}

		/// <summary>
		/// Random variable with a gaussian distribution.
		/// </summary>
		Vector2 GaussianRandomVariable() 
		{
			float x1, x2, w;
			do 
			{
				x1 = 2.0f * UnityEngine.Random.value - 1.0f;
				x2 = 2.0f * UnityEngine.Random.value - 1.0f;
				w = x1 * x1 + x2 * x2;
			} 
			while ( w >= 1.0f );
			
			w = Mathf.Sqrt((-2.0f * Mathf.Log(w)) / w);
			return new Vector2(x1 * w, x2 * w);
		}

		/// <summary>
		/// Gets the spectrum vaule for grid position n,m.
		/// </summary>
		float PhillipsSpectrum(int n_prime, int m_prime) 
		{
			Vector2 k = new Vector2(Mathf.PI * (2 * n_prime - N) / m_length, Mathf.PI * (2 * m_prime - N) / m_length);
			float k_length  = k.magnitude;
			if (k_length < 0.000001f) return 0.0f;
			
			float k_length2 = k_length  * k_length;
			float k_length4 = k_length2 * k_length2;
			
			k.Normalize();
			
			float k_dot_w   = Vector2.Dot(k, m_windDirection);
			float k_dot_w2  = k_dot_w * k_dot_w * k_dot_w * k_dot_w * k_dot_w * k_dot_w;
			
			float w_length  = m_windSpeed.magnitude;
			float L         = w_length * w_length / GRAVITY;
			float L2        = L * L;
			
			float damping   = 0.001f;
			float l2        = L2 * damping * damping;
			
			return m_waveAmp * Mathf.Exp(-1.0f / (k_length2 * L2)) / k_length4 * k_dot_w2 * Mathf.Exp(-k_length2 * l2);
		}
		
		float Dispersion(int n_prime, int m_prime) 
		{
			float w_0 = 2.0f * Mathf.PI / 200.0f;
			float kx = Mathf.PI * (2 * n_prime - N) / m_length;
			float kz = Mathf.PI * (2 * m_prime - N) / m_length;
			return Mathf.Floor(Mathf.Sqrt(GRAVITY * Mathf.Sqrt(kx * kx + kz * kz)) / w_0) * w_0;
		}

		/// <summary>
		/// Inits the spectrum for time period t.
		/// </summary>
		Vector2 InitSpectrum(float t, int n_prime, int m_prime) 
		{
			int index = m_prime * Nplus1 + n_prime;
		
			float omegat = m_dispersionTable[index] * t;
		
			float cos = Mathf.Cos(omegat);
			float sin = Mathf.Sin(omegat);
		
			float c0a = m_spectrum[index].x*cos - m_spectrum[index].y*sin;
			float c0b = m_spectrum[index].x*sin + m_spectrum[index].y*cos;
			
			float c1a = m_spectrum_conj[index].x*cos - m_spectrum_conj[index].y*-sin;
			float c1b = m_spectrum_conj[index].x*-sin + m_spectrum_conj[index].y*cos;
			
			return new Vector2(c0a+c1a, c0b+c1b);
		}

		/// <summary>
		/// Runs a threaded task
		/// </summary>
		void RunThreaded(object o)
		{
			
			Nullable<float> time  = o as Nullable<float>;

			EvaluateWavesFFT(time.Value);

			done = true;
			
		}

		/// <summary>
		/// Evaluates the waves for time period t. Must be thread safe.
		/// </summary>
		void EvaluateWavesFFT(float t) 
		{


			float kx, kz, len, lambda = -1.0f;
			int index, index1;
		
			for (int m_prime = 0; m_prime < N; m_prime++) 
			{
				kz = Mathf.PI * (2.0f * m_prime - N) / m_length;
				
				for (int n_prime = 0; n_prime < N; n_prime++) 
				{
					kx = Mathf.PI*(2 * n_prime - N) / m_length;
					len = Mathf.Sqrt(kx * kx + kz * kz);
					index = m_prime * N + n_prime;
					
					Vector2 c = InitSpectrum(t, n_prime, m_prime);
		
					m_heightBuffer[1,index].x = c.x;
					m_heightBuffer[1,index].y = c.y;
					
					m_slopeBuffer[1,index].x = -c.y*kx;
					m_slopeBuffer[1,index].y = c.x*kx;
					
					m_slopeBuffer[1,index].z = -c.y*kz;
					m_slopeBuffer[1,index].w = c.x*kz;
					
					if (len < 0.000001f) 
					{
						m_displacementBuffer[1,index].x = 0.0f;
						m_displacementBuffer[1,index].y = 0.0f;
						m_displacementBuffer[1,index].z = 0.0f;
						m_displacementBuffer[1,index].w = 0.0f;
					} 
					else 
					{
						m_displacementBuffer[1,index].x = -c.y * -(kx/len);
						m_displacementBuffer[1,index].y = c.x * -(kx/len);
						m_displacementBuffer[1,index].z = -c.y * -(kz/len);
						m_displacementBuffer[1,index].w = c.x * -(kz/len);
					}	
				}
			}
			
			m_fourier.PeformFFT(0, m_heightBuffer, m_slopeBuffer, m_displacementBuffer);
	
			int sign;
			float[] signs = new float[]{ 1.0f, -1.0f };
			Vector3 n;
			
			for (int m_prime = 0; m_prime < N; m_prime++) 
			{
				for (int n_prime = 0; n_prime < N; n_prime++) 
				{
					index  = m_prime * N + n_prime;			// index into buffers
					index1 = m_prime * Nplus1 + n_prime;	// index into vertices
		
					sign = (int)signs[(n_prime + m_prime) & 1];
		
					// height
					m_vertices[index1].y = m_heightBuffer[1, index].x * sign;
		
					// displacement
					m_vertices[index1].x = m_position[index1].x + m_displacementBuffer[1, index].x * lambda * sign;
					m_vertices[index1].z = m_position[index1].z + m_displacementBuffer[1, index].z * lambda * sign;
					
					// normal
					n = new Vector3(-m_slopeBuffer[1, index].x * sign, 1.0f, -m_slopeBuffer[1, index].z * sign);
					n.Normalize();
					
					m_normals[index1].x =  n.x;
					m_normals[index1].y =  n.y;
					m_normals[index1].z =  n.z;
		
					// for tiling
					if (n_prime == 0 && m_prime == 0) 
					{
						m_vertices[index1 + N + Nplus1 * N].y = m_heightBuffer[1, index].x * sign;
		
						m_vertices[index1 + N + Nplus1 * N].x = m_position[index1 + N + Nplus1 * N].x + m_displacementBuffer[1, index].x * lambda * sign;
						m_vertices[index1 + N + Nplus1 * N].z = m_position[index1 + N + Nplus1 * N].z + m_displacementBuffer[1, index].z * lambda * sign;
					
						m_normals[index1 + N + Nplus1 * N].x =  n.x;
						m_normals[index1 + N + Nplus1 * N].y =  n.y;
						m_normals[index1 + N + Nplus1 * N].z =  n.z;
					}
					if (n_prime == 0) 
					{
						m_vertices[index1 + N].y = m_heightBuffer[1, index].x * sign;
		
						m_vertices[index1 + N].x = m_position[index1 + N].x + m_displacementBuffer[1, index].x * lambda * sign;
						m_vertices[index1 + N].z = m_position[index1 + N].z + m_displacementBuffer[1, index].z * lambda * sign;
					
						m_normals[index1 + N].x =  n.x;
						m_normals[index1 + N].y =  n.y;
						m_normals[index1 + N].z =  n.z;
					}
					if (m_prime == 0) 
					{
						m_vertices[index1 + Nplus1 * N].y = m_heightBuffer[1, index].x * sign;
		
						m_vertices[index1 + Nplus1 * N].x = m_position[index1 + Nplus1 * N].x + m_displacementBuffer[1, index].x * lambda * sign;
						m_vertices[index1 + Nplus1 * N].z = m_position[index1 + Nplus1 * N].z + m_displacementBuffer[1, index].z * lambda * sign;
					
						m_normals[index1 + Nplus1 * N].x =  n.x;
						m_normals[index1 + Nplus1 * N].y =  n.y;
						m_normals[index1 + Nplus1 * N].z =  n.z;
					}
				}
			}



		}

		Mesh MakeMesh(int size) 
		{
		
			Vector3[] vertices = new Vector3[size*size];
			Vector2[] texcoords = new Vector2[size*size];
			Vector3[] normals = new Vector3[size*size];
			int[] indices = new int[size*size*6];
			
			for(int x = 0; x < size; x++)
			{
				for(int y = 0; y < size; y++)
				{
					Vector2 uv = new Vector3( (float)x / (float)(size-1), (float)y / (float)(size-1) );
					Vector3 pos = new Vector3(x, 0.0f, y);
					Vector3 norm = new Vector3(0.0f, 1.0f, 0.0f);
					
					texcoords[x+y*size] = uv;
					vertices[x+y*size] = pos;
					normals[x+y*size] = norm;
				}
			}
			
			int num = 0;
			for(int x = 0; x < size-1; x++)
			{
				for(int y = 0; y < size-1; y++)
				{
					indices[num++] = x + y * size;
					indices[num++] = x + (y+1) * size;
					indices[num++] = (x+1) + y * size;
			
					indices[num++] = x + (y+1) * size;
					indices[num++] = (x+1) + (y+1) * size;
					indices[num++] = (x+1) + y * size;
				}
			}
			
			Mesh mesh = new Mesh();
		
			mesh.vertices = vertices;
			mesh.uv = texcoords;
			mesh.triangles = indices;
			mesh.normals = normals;
			
			return mesh;
		}
	}

}

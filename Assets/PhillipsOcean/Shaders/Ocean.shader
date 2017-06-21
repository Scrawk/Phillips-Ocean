Shader "Ocean/Ocean" 
{
	Properties 
	{
		_SeaColor("SeaColor", Color) = (1,1,1,1)
		_SkyColor("SkyColor", Color) = (1,1,1,1)
	}
	SubShader 
	{
		Tags { "RenderType"="Opaque" }
		LOD 200
		
		CGPROGRAM
		#pragma surface surf Lambert

		float3 _SeaColor, _SkyColor;
		sampler2D _FresnelLookUp;
		
		struct Input 
		{
			float3 worldNormal;
			float3 worldPos;
			float3 worldRefl;
		};
		
		float Fresnel(float3 V, float3 N)
		{
			float costhetai = abs(dot(V, N));
			return tex2D(_FresnelLookUp, float2(costhetai, 0.0)).a;
		}

		void surf (Input IN, inout SurfaceOutput o) 
		{
			float3 V = normalize(_WorldSpaceCameraPos-IN.worldPos);
			float3 N = IN.worldNormal;
		
			float fresnel = Fresnel(V, N);
			
			o.Albedo = lerp(_SeaColor, _SkyColor, fresnel);
			o.Alpha = 1.0;
		}
		ENDCG
	} 
	FallBack "Diffuse"
}
















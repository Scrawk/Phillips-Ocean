# Phillips-Ocean

This is a ocean project  using Phillips spectrum to generate the waves. Some time around 2001 Jerry Tessendorf released a paper called 'Simulating Ocean Water' and in it he out lined all the methods and math need for simulating oceans. In this paper the math for using Phillips spectrum was outlined and I came across a [blog](http://www.keithlantz.net/2011/11/ocean-simulation-part-two-using-the-fast-fourier-transform/) that had converted the math to code.

His project was written in C++ so it was just a matter of converting that to a C# script for Unity. I have left all the math the same but have made some changes and restructured the code a bit and the Fourier transform now runs on its own thread resulting in a big performance increase.

You can download a Unity package [here](https://app.box.com/s/mlqw507k7nz7vvmbn9ttjvxo7e5981to).

![Phillips Ocean](https://static.wixstatic.com/media/1e04d5_07fba343a78f48c996b353a170c05be5~mv2.png/v1/fill/w_486,h_486,al_c,usm_0.66_1.00_0.01/1e04d5_07fba343a78f48c996b353a170c05be5~mv2.png)

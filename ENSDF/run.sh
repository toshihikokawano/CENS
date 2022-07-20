#!/bin/sh

src=orig
dst=adopted

for file in\
 ensdf.001 ensdf.002 ensdf.003 ensdf.004 ensdf.005 ensdf.006 ensdf.007 ensdf.008 ensdf.009 ensdf.010\
 ensdf.011 ensdf.012 ensdf.013 ensdf.014 ensdf.015 ensdf.016 ensdf.017 ensdf.018 ensdf.019 ensdf.020\
 ensdf.021 ensdf.022 ensdf.023 ensdf.024 ensdf.025 ensdf.026 ensdf.027 ensdf.028 ensdf.029 ensdf.030\
 ensdf.031 ensdf.032 ensdf.033 ensdf.034 ensdf.035 ensdf.036 ensdf.037 ensdf.038 ensdf.039 ensdf.040\
 ensdf.041 ensdf.042 ensdf.043 ensdf.044 ensdf.045 ensdf.046 ensdf.047 ensdf.048 ensdf.049 ensdf.050\
 ensdf.051 ensdf.052 ensdf.053 ensdf.054 ensdf.055 ensdf.056 ensdf.057 ensdf.058 ensdf.059 ensdf.060\
 ensdf.061 ensdf.062 ensdf.063 ensdf.064 ensdf.065 ensdf.066 ensdf.067 ensdf.068 ensdf.069 ensdf.070\
 ensdf.071 ensdf.072 ensdf.073 ensdf.074 ensdf.075 ensdf.076 ensdf.077 ensdf.078 ensdf.079 ensdf.080\
 ensdf.081 ensdf.082 ensdf.083 ensdf.084 ensdf.085 ensdf.086 ensdf.087 ensdf.088 ensdf.089 ensdf.090\
 ensdf.091 ensdf.092 ensdf.093 ensdf.094 ensdf.095 ensdf.096 ensdf.097 ensdf.098 ensdf.099 ensdf.100\
 ensdf.101 ensdf.102 ensdf.103 ensdf.104 ensdf.105 ensdf.106 ensdf.107 ensdf.108 ensdf.109 ensdf.110\
 ensdf.111 ensdf.112 ensdf.113 ensdf.114 ensdf.115 ensdf.116 ensdf.117 ensdf.118 ensdf.119 ensdf.120\
 ensdf.121 ensdf.122 ensdf.123 ensdf.124 ensdf.125 ensdf.126 ensdf.127 ensdf.128 ensdf.129 ensdf.130\
 ensdf.131 ensdf.132 ensdf.133 ensdf.134 ensdf.135 ensdf.136 ensdf.137 ensdf.138 ensdf.139 ensdf.140\
 ensdf.141 ensdf.142 ensdf.143 ensdf.144 ensdf.145 ensdf.146 ensdf.147 ensdf.148 ensdf.149 ensdf.150\
 ensdf.151 ensdf.152 ensdf.153 ensdf.154 ensdf.155 ensdf.156 ensdf.157 ensdf.158 ensdf.159 ensdf.160\
 ensdf.161 ensdf.162 ensdf.163 ensdf.164 ensdf.165 ensdf.166 ensdf.167 ensdf.168 ensdf.169 ensdf.170\
 ensdf.171 ensdf.172 ensdf.173 ensdf.174 ensdf.175 ensdf.176 ensdf.177 ensdf.178 ensdf.179 ensdf.180\
 ensdf.181 ensdf.182 ensdf.183 ensdf.184 ensdf.185 ensdf.186 ensdf.187 ensdf.188 ensdf.189 ensdf.190\
 ensdf.191 ensdf.192 ensdf.193 ensdf.194 ensdf.195 ensdf.196 ensdf.197 ensdf.198 ensdf.199 ensdf.200\
 ensdf.201 ensdf.202 ensdf.203 ensdf.204 ensdf.205 ensdf.206 ensdf.207 ensdf.208 ensdf.209 ensdf.210\
 ensdf.211 ensdf.212 ensdf.213 ensdf.214 ensdf.215 ensdf.216 ensdf.217 ensdf.218 ensdf.219 ensdf.220\
 ensdf.221 ensdf.222 ensdf.223 ensdf.224 ensdf.225 ensdf.226 ensdf.227 ensdf.228 ensdf.229 ensdf.230\
 ensdf.231 ensdf.232 ensdf.233 ensdf.234 ensdf.235 ensdf.236 ensdf.237 ensdf.238 ensdf.239 ensdf.240\
 ensdf.241 ensdf.242 ensdf.243 ensdf.244 ensdf.245 ensdf.246 ensdf.247 ensdf.248 ensdf.249 ensdf.250\
 ensdf.251 ensdf.252 ensdf.253 ensdf.254 ensdf.255 ensdf.256 ensdf.257 ensdf.258 ensdf.259 ensdf.260\
 ensdf.261 ensdf.262 ensdf.263 ensdf.264 ensdf.265 ensdf.266 ensdf.267 ensdf.268 ensdf.269 ensdf.270\
 ensdf.271 ensdf.272 ensdf.273 ensdf.274 ensdf.275 ensdf.276 ensdf.277 ensdf.278 ensdf.279 ensdf.280\
 ensdf.281 ensdf.282 ensdf.283 ensdf.284 ensdf.285 ensdf.286 ensdf.287 ensdf.288 ensdf.289 ensdf.290\
 ensdf.291 ensdf.292 ensdf.293 ensdf.294 ensdf.295 ensdf.296 ensdf.297 ensdf.298 ensdf.299 ensdf.300; do
 name=$src/$file
 echo $name
 ./adoptedlevel.pl $name
 mv ENSDF*.dat $dst
done


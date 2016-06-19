# NolinearTimeSeriesAnalysis
USAGE
The codes in the toolbox can be used to perform nonlinear time series analysis on single(or multi) channel data. This is done by mapping the single channel data to phase space representation using Taken's embedding theorem (compute_psv.m). The parameters - optimal delay and dimension are estimated using first minimum of MI (compute_tau.m) and FNN method (compute_dim) respectively. The recurrence network can be constructed from the phase space vector using ComputeRecurrenceNetwork_ANN.m or ComputeRecurrenceNetwork_fixedRR.m. The topology of the RN can be further analysed using graph theoreticl quantifiers (you need BCT toolbox for this). One can also compute the complexity-entrropy plane using get_mpr_complexity.m for which the ordinal patterns are computed using get_ordinal_pattern_dist.m (see the function descp for more details). Also, the tool box contains python codes to generate variety of uni(or multi) variate surrogate data.


LICENSE INFORMATION:

NonlinearTimeSeriesAnalysis toolbox is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

NonlinearTimeSeriesAnalysis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

Narayan P Subramaniyam - npsubramaniyam@gmail.com (01.11.2015)

The following publication must be cited in case the codes from the toolbox are used : 

Signatures of chaotic and stochastic dynamics uncovered with ε-recurrence networks, NP Subramaniyam, JF Donges, J Hyttinen - Proc. R. Soc. A, 2015

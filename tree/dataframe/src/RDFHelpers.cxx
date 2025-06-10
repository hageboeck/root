// Author: Stefan Wunsch, Enrico Guiraud CERN  09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDFHelpers.hxx"
#include "TROOT.h"      // IsImplicitMTEnabled
#include "TError.h"     // Warning
#include "TStopwatch.h"
#include "RConfigure.h" // R__USE_IMT
#include "ROOT/RLogger.hxx"
#include "ROOT/RDF/RLoopManager.hxx" // for RLoopManager
#include "ROOT/RDF/Utils.hxx"
#include "ROOT/RResultHandle.hxx"    // for RResultHandle, RunGraphs
#include "ROOT/RSlotStack.hxx"
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif // R__USE_IMT

#include <TCanvas.h>
#include <TPRegexp.h>

#include <algorithm>
#include <future>
#include <iostream>
#include <set>
#include <cstdio>

// TODO, this function should be part of core libraries
#include <numeric>
#if (!defined(_WIN32)) && (!defined(_WIN64))
#include <unistd.h>
#endif

#if defined(_WIN32) || defined(_WIN64)
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <io.h>
#include <Windows.h>
#else
#include <sys/ioctl.h>
#endif

// Get terminal size for progress bar
int get_tty_size()
{
#if defined(_WIN32) || defined(_WIN64)
   if (!_isatty(_fileno(stdout)))
      return 0;
   int width = 0;
   CONSOLE_SCREEN_BUFFER_INFO csbi;
   if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
      width = (int)(csbi.srWindow.Right - csbi.srWindow.Left + 1);
   return width;
#else
   int width = 0;
   struct winsize w;
   ioctl(fileno(stdout), TIOCGWINSZ, &w);
   width = (int)(w.ws_col);
   return width;
#endif
}

using ROOT::RDF::RResultHandle;

unsigned int ROOT::RDF::RunGraphs(std::vector<RResultHandle> handles)
{
   if (handles.empty()) {
      Warning("RunGraphs", "Got an empty list of handles, now quitting.");
      return 0u;
   }

   // Check that there are results which have not yet been run
   const unsigned int nToRun =
      std::count_if(handles.begin(), handles.end(), [](const auto &h) { return !h.IsReady(); });
   if (nToRun < handles.size()) {
      Warning("RunGraphs", "Got %zu handles from which %zu link to results which are already ready.", handles.size(),
              handles.size() - nToRun);
   }
   if (nToRun == 0u)
      return 0u;

   // Find the unique event loops
   auto sameGraph = [](const RResultHandle &a, const RResultHandle &b) { return a.fLoopManager < b.fLoopManager; };
   std::set<RResultHandle, decltype(sameGraph)> s(handles.begin(), handles.end(), sameGraph);
   std::vector<RResultHandle> uniqueLoops(s.begin(), s.end());

   // Trigger jitting. One call is enough to jit the code required by all computation graphs.
   TStopwatch sw;
   sw.Start();
   {
      const auto effectiveVerbosity = ROOT::Internal::GetChannelOrManager(ROOT::Detail::RDF::RDFLogChannel())
                                         .GetEffectiveVerbosity(ROOT::RLogManager::Get());
      if (effectiveVerbosity >= ROOT::ELogLevel::kDebug + 10) {
         // a very high verbosity was requested, let's not silence anything
         uniqueLoops[0].fLoopManager->Jit();
      } else {
         // silence logs from RLoopManager::Jit: RunGraphs does its own logging
         auto silenceRDFLogs = ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::ELogLevel::kError);
         uniqueLoops[0].fLoopManager->Jit();
      }
   }
   sw.Stop();
   R__LOG_INFO(ROOT::Detail::RDF::RDFLogChannel())
      << "Just-in-time compilation phase for RunGraphs (" << uniqueLoops.size()
      << " unique computation graphs) completed"
      << (sw.RealTime() > 1e-3 ? " in " + std::to_string(sw.RealTime()) + " seconds." : " in less than 1ms.");

   // Trigger the unique event loops
   auto slotStack = std::make_shared<ROOT::Internal::RSlotStack>(ROOT::GetThreadPoolSize());
   auto run = [&slotStack](RResultHandle &h) {
      if (h.fLoopManager) {
         h.fLoopManager->SetSlotStack(slotStack);
         h.fLoopManager->Run(/*jit=*/false);
      }
   };

   sw.Start();
#ifdef R__USE_IMT
   if (ROOT::IsImplicitMTEnabled()) {
      ROOT::TThreadExecutor{}.Foreach(run, uniqueLoops);
   } else {
#endif
      std::for_each(uniqueLoops.begin(), uniqueLoops.end(), run);
#ifdef R__USE_IMT
   }
#endif
   sw.Stop();
   R__LOG_INFO(ROOT::Detail::RDF::RDFLogChannel())
      << "Finished RunGraphs run (" << uniqueLoops.size() << " unique computation graphs, " << sw.CpuTime() << "s CPU, "
      << sw.RealTime() << "s elapsed).";

   return uniqueLoops.size();
}

/// Draw an interactive overview of a dataset using RDataFrame.
/// \param treename Name of a tree or RNTuple containing the dataset.
/// \param filenameglob Glob expression specifying the files of this dataset.
/// \param columns Regular expression selecting the columns of the dataset to be drawn. Leaving empty draws all columns.
/// \param selection Cuts to be applied when drawing the overview.
/// \param events TODO: Implement event range
void ROOT::RDF::Draw(std::string treename, std::string filenameglob, std::string columns, std::string selection, ULong64_t events) {
   using namespace std::chrono_literals;

   struct BranchData {
      std::string name;
      double min;
      double max;
      std::unique_ptr<TH1D> partial;
      ROOT::RDF::RResultPtr<TH1D> resultPtr;
      std::unique_ptr<std::mutex> mutex{new std::mutex};
   };
   std::vector<BranchData> histos;
   static std::vector<std::shared_ptr<TH1D>> histoLifeline;
   histoLifeline.clear();


   {
      // TODO: Don't use range, and use MT
      auto rangedRdf = ROOT::RDataFrame{treename, filenameglob}.Range(10000);
      auto colNameTmp = rangedRdf.GetColumnNames();
      std::unique_ptr<TPRegexp> regex;
      if (!columns.empty()) {
        regex = std::make_unique<TPRegexp>(columns.c_str());
      }
      std::ostringstream buf;
      buf << std::left;
      const auto rows = (colNameTmp.size()+3) / 4;
      std::vector<std::string> columnNames;

      for (unsigned int row = 0; row < rows; ++row) {
         for (auto index = row; index < colNameTmp.size(); index += rows) {
            auto const & colName = colNameTmp[index];
            bool match = regex ? regex->MatchB(colName.c_str()) : true;
            if (match) columnNames.push_back(colName);
            buf << (match ? "[x]" : "[ ]") << std::setw(30) << colName;
         }
         buf << "\n";
      }
      std::cout << buf.str();

      std::vector<std::pair<ROOT::RDF::RResultPtr<double>,ROOT::RDF::RResultPtr<double>>> minMaxResPtrs;
      minMaxResPtrs.reserve(columnNames.size());
      for (auto const & colname : columnNames) {
         minMaxResPtrs.emplace_back(rangedRdf.Min(colname), rangedRdf.Max(colname));
      }
      histos.reserve(columnNames.size());
      for (unsigned int i = 0; i < columnNames.size(); ++i) {
         const auto min = minMaxResPtrs[i].first.GetValue(), max = minMaxResPtrs[i].second.GetValue();
         histos.emplace_back(columnNames[i], min, 1.1*max);
      }
   }

   std::sort(histos.begin(), histos.end(), [](BranchData const & lhs, BranchData const & rhs){ return lhs.name < rhs.name; });

   ROOT::EnableImplicitMT();
   auto future = std::async(std::launch::async, [&]() -> void {
      ROOT::RDataFrame rdf(treename, filenameglob);
      ROOT::RDF::RNode rootNode{rdf};
      if (!selection.empty()) {
        rootNode = rdf.Filter(selection);
      }

      for (unsigned int i = 0; i < histos.size(); ++i) {
         BranchData & hw = histos[i];
         const bool validLimits = hw.min < hw.max;
         const std::string col = hw.name;
         hw.resultPtr = rootNode.Histo1D(TH1DModel{col.c_str(), (col + ";" + col + ";Events").c_str(), 100,
               validLimits ? hw.min : 0., validLimits ? hw.max : 1.1}, col);
         hw.resultPtr.OnPartialResult(50000, [&hw](TH1D const & partial){
            std::scoped_lock lock{*hw.mutex};
            if (!hw.partial) {
               hw.partial.reset(new TH1D{partial});
            } else {
               hw.partial->Add(&partial);
            }
         });
      }

      histos.front().resultPtr.GetPtr();
   });


   TCanvas * c = new TCanvas{("RDF: " + treename).c_str(), ("RDataFrame overview of " + treename).c_str()};
   c->SetCanvasSize(2048, 2048);
   c->SetWindowSize(2048, 1024);

   const unsigned int cols = std::sqrt(histos.size());
   c->Divide(cols, histos.size()/cols);

   while (future.wait_for(2s) == std::future_status::timeout) {
      for (unsigned int i = 0; i < histos.size(); ++i) {
         auto pad = c->cd(i+1);
         BranchData & hw = histos[i];
         std::scoped_lock lock{*hw.mutex};
         if (hw.partial) {
            hw.partial->DrawCopy();
            pad->Modified();
            pad->Update();
         }
      }
      c->Update();
   }

   future.get();
   for (unsigned int i = 0; i < histos.size(); ++i) {
      c->cd(i+1)->Clear();
      BranchData & hw = histos[i];
      hw.partial = nullptr;
      hw.resultPtr->Draw();
      histoLifeline.push_back(hw.resultPtr.GetSharedPtr());
   }
   c->Draw();
}

namespace ROOT::RDF::Experimental {

/// Add systematic variations to a snapshot.
/// \param[in] resPtr The snapshot instance for which variations should be produced.
/// \return A \ref ROOT::RDF::RResultPtr to a new dataframe that has nominal and varied columns.
///
/// VariationsFor does not trigger the event loop. The event loop is only triggered
/// upon first access to a valid key, similarly to what happens with RResultPtr.
///
/// See RDataFrame's \ref ROOT::RDF::RInterface::Vary() "Vary" method for more information and example usages.
RResultPtr<SnapshotResult_t> VariationsFor(RResultPtr<SnapshotResult_t> resPtr)
{
   R__ASSERT(resPtr != nullptr && "Calling VariationsFor on an empty RResultPtr");

   // populate parts of the computation graph for which we only have "empty shells", e.g. RJittedActions and
   // RJittedFilters
   resPtr.fLoopManager->Jit();

   std::unique_ptr<RDFInternal::RActionBase> variedAction;

   std::shared_ptr<RDFInternal::RActionBase> nominalAction = resPtr.fActionPtr;
   std::vector<std::string> variations = nominalAction->GetVariations();
   const auto nVariations = variations.size();

   if (nVariations > 0) {
      // Create the RVariedAction and inject it in the computation graph.
      // This recursively creates all the required varied column readers and upstream nodes of the computation graph.
      variedAction = nominalAction->MakeVariedAction(std::vector<void*>(nVariations));
   }

   return ROOT::Detail::RDF::MakeResultPtr<SnapshotResult_t>(resPtr.fObjPtr, *resPtr.fLoopManager, std::move(variedAction));
}

void ThreadsPerTH3(unsigned int N)
{
   ROOT::Internal::RDF::NThreadPerTH3() = N;
}

ProgressHelper::ProgressHelper(std::size_t increment, unsigned int totalFiles, unsigned int progressBarWidth,
                               unsigned int printInterval, bool useColors)
   : fPrintInterval(printInterval),
     fIncrement{increment},
     fBarWidth{progressBarWidth = int(get_tty_size() / 4)},
     fTotalFiles{totalFiles},
#if defined(_WIN32) || defined(_WIN64)
     fIsTTY{_isatty(_fileno(stdout)) != 0},
     fUseShellColours{false && useColors}
#else
     fIsTTY{isatty(fileno(stdout)) == 1},
     fUseShellColours{useColors && fIsTTY} // Control characters only with terminals.
#endif
{
}

/// Compute a running mean of events/s.
double ProgressHelper::EvtPerSec() const
{
   if (fEventsPerSecondStatisticsIndex < fEventsPerSecondStatistics.size())
      return std::accumulate(fEventsPerSecondStatistics.begin(),
                             fEventsPerSecondStatistics.begin() + fEventsPerSecondStatisticsIndex, 0.) /
             fEventsPerSecondStatisticsIndex;
   else
      return std::accumulate(fEventsPerSecondStatistics.begin(), fEventsPerSecondStatistics.end(), 0.) /
             fEventsPerSecondStatistics.size();
}

/// Record current event counts and time stamp, populate evts/s statistics array.
std::pair<std::size_t, std::chrono::seconds> ProgressHelper::RecordEvtCountAndTime()
{
   using namespace std::chrono;

   auto currentEventCount = fProcessedEvents.load();
   auto eventsPerTimeInterval = currentEventCount - fLastProcessedEvents;
   fLastProcessedEvents = currentEventCount;

   auto oldPrintTime = fLastPrintTime;
   auto newPrintTime = system_clock::now();
   fLastPrintTime = newPrintTime;

   duration<double> secondsCurrentInterval = newPrintTime - oldPrintTime;
   fEventsPerSecondStatistics[fEventsPerSecondStatisticsIndex++ % fEventsPerSecondStatistics.size()] =
      eventsPerTimeInterval / secondsCurrentInterval.count();

   return {currentEventCount, duration_cast<seconds>(newPrintTime - fBeginTime)};
}

namespace {

struct RestoreStreamState {
   RestoreStreamState(std::ostream &stream) : fStream(stream), fFlags(stream.flags()), fFillChar(stream.fill()) {}
   ~RestoreStreamState()
   {
      fStream.flags(fFlags);
      fStream.fill(fFillChar);
   }

   std::ostream &fStream;
   std::ios_base::fmtflags fFlags;
   std::ostream::char_type fFillChar;
};

/// Format std::chrono::seconds as `1:30m`.
std::ostream &operator<<(std::ostream &stream, std::chrono::seconds elapsedSeconds)
{
   RestoreStreamState restore(stream);
   auto h = std::chrono::duration_cast<std::chrono::hours>(elapsedSeconds);
   auto m = std::chrono::duration_cast<std::chrono::minutes>(elapsedSeconds - h);
   auto s = (elapsedSeconds - h - m).count();

   if (h.count() > 0)
      stream << h.count() << ':' << std::setw(2) << std::right << std::setfill('0');
   stream << m.count() << ':' << std::setw(2) << std::right << std::setfill('0') << s;
   return stream << (h.count() > 0 ? 'h' : 'm');
}

} // namespace

/// Print event and time statistics.
void ProgressHelper::PrintStats(std::ostream &stream, std::size_t currentEventCount,
                                std::chrono::seconds elapsedSeconds) const
{
   RestoreStreamState restore(stream);
   auto evtpersec = EvtPerSec();
   auto GetNEventsOfCurrentFile = ComputeNEventsSoFar();
   auto currentFileIdx = ComputeCurrentFileIdx();
   auto totalFiles = fTotalFiles;

   if (fUseShellColours)
      stream << "\033[35m";
   stream << "["
          << "Elapsed time: " << elapsedSeconds << "  ";
   if (fUseShellColours)
      stream << "\033[0m";
   stream << "processing file: " << currentFileIdx << " / " << totalFiles << "  ";

   // Event counts:
   if (fUseShellColours)
      stream << "\033[32m";

   stream << "processed evts: " << currentEventCount;
   if (GetNEventsOfCurrentFile != 0) {
      stream << " / " << std::scientific << std::setprecision(2) << GetNEventsOfCurrentFile;
   }
   stream << "  ";

   if (fUseShellColours)
      stream << "\033[0m";

   // events/s
   stream << std::scientific << std::setprecision(2) << evtpersec << " evt/s";

   // Time statistics:
   if (GetNEventsOfCurrentFile != 0) {
      if (fUseShellColours)
         stream << "\033[35m";
      std::chrono::seconds remainingSeconds(
         static_cast<long long>((ComputeNEventsSoFar() - currentEventCount) / evtpersec));
      stream << " " << remainingSeconds << " "
             << " remaining time (per file being processed)";
      if (fUseShellColours)
         stream << "\033[0m";
   }

   stream << "]   ";
}

void ProgressHelper::PrintStatsFinal(std::ostream &stream, std::chrono::seconds elapsedSeconds) const
{
   RestoreStreamState restore(stream);
   auto totalEvents = ComputeNEventsSoFar();
   auto totalFiles = fTotalFiles;

   if (fUseShellColours)
      stream << "\033[35m";
   stream << "["
          << "Total elapsed time: " << elapsedSeconds << "  ";
   if (fUseShellColours)
      stream << "\033[0m";
   stream << "processed files: " << totalFiles << " / " << totalFiles << "  ";

   // Event counts:
   if (fUseShellColours)
      stream << "\033[32m";

   stream << "processed evts: " << totalEvents;
   if (totalEvents != 0) {
      stream << " / " << std::scientific << std::setprecision(2) << totalEvents;
   }

   if (fUseShellColours)
      stream << "\033[0m";

   stream << "]   ";
}

/// Print a progress bar of width `ProgressHelper::fBarWidth` if `fGetNEventsOfCurrentFile` is known.
void ProgressHelper::PrintProgressBar(std::ostream &stream, std::size_t currentEventCount) const
{
   auto GetNEventsOfCurrentFile = ComputeNEventsSoFar();
   if (GetNEventsOfCurrentFile == 0)
      return;

   RestoreStreamState restore(stream);

   double completion = double(currentEventCount) / GetNEventsOfCurrentFile;
   unsigned int nBar = std::min(completion, 1.) * fBarWidth;

   std::string bars(std::max(nBar, 1u), '=');
   bars.back() = (nBar == fBarWidth) ? '=' : '>';

   if (fUseShellColours)
      stream << "\033[33m";
   stream << '|' << std::setfill(' ') << std::setw(fBarWidth) << std::left << bars << "|   ";
   if (fUseShellColours)
      stream << "\033[0m";
}
//*/

class ProgressBarAction final : public ROOT::Detail::RDF::RActionImpl<ProgressBarAction> {
public:
   using Result_t = int;

private:
   std::shared_ptr<ProgressHelper> fHelper;
   std::shared_ptr<int> fDummyResult = std::make_shared<int>();

public:
   ProgressBarAction(std::shared_ptr<ProgressHelper> r) : fHelper(std::move(r)) {}

   std::shared_ptr<Result_t> GetResultPtr() const { return fDummyResult; }

   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}

   void Exec(unsigned int) {}

   void Finalize()
   {
      std::mutex fPrintMutex;
      if (!fPrintMutex.try_lock())
         return;
      std::lock_guard<std::mutex> lockGuard(fPrintMutex, std::adopt_lock);
      const auto &[eventCount, elapsedSeconds] = fHelper->RecordEvtCountAndTime();

      // The next line resets the current line output in the terminal.
      // Brings the cursor at the beginning ('\r'), prints whitespace with the
      // same length as the terminal size, then resets the cursor again so we
      // can print the final stats on a clean line.
      std::cout << '\r' << std::string(get_tty_size(), ' ') << '\r';
      fHelper->PrintStatsFinal(std::cout, elapsedSeconds);
      std::cout << '\n';
   }

   std::string GetActionName() { return "ProgressBar"; }
   // dummy implementation of PartialUpdate
   int &PartialUpdate(unsigned int) { return *fDummyResult; }

   ROOT::RDF::SampleCallback_t GetSampleCallback() final
   {
      return [this](unsigned int slot, const ROOT::RDF::RSampleInfo &id) {
         this->fHelper->registerNewSample(slot, id);
         return this->fHelper->ComputeNEventsSoFar();
      };
   }
};

void AddProgressBar(ROOT::RDF::RNode node)
{
   auto total_files = node.GetNFiles();
   auto progress = std::make_shared<ProgressHelper>(1000, total_files);
   ProgressBarAction c(progress);
   auto r = node.Book<>(c);
   r.OnPartialResultSlot(1000, [progress](unsigned int slot, auto &&arg) { (*progress)(slot, arg); });
}

void AddProgressBar(ROOT::RDataFrame dataframe)
{
   auto node = ROOT::RDF::AsRNode(dataframe);
   ROOT::RDF::Experimental::AddProgressBar(node);
}

} // namespace ROOT

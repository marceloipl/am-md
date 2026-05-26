# GTLab Toolbox - Documentação Final Consolidada (PT)

Este pacote contém a versão portuguesa do relatório final da GTLab Toolbox.

## Conteúdo

- `final_report.tex` - documento LaTeX mestre em português.
- `final_report.pdf` - versão PDF compilada do relatório final.
- `chapters/analysis/` - capítulos de análise e documentação dos ficheiros MATLAB/Octave.
- `chapters/testing/controlled_testing.tex` - testes controlados, resultados e conclusões.
- `src/GTLabtools/` - código-fonte MATLAB/Octave analisado.
- `original_documents/` - documentos originais usados como base, preservados para rastreabilidade.

## Como compilar

```bash
pdflatex final_report.tex
pdflatex final_report.tex
```

## Como reproduzir os testes

1. Abrir MATLAB ou GNU Octave.
2. Adicionar a pasta de código: `addpath('src/GTLabtools')`.
3. Executar os blocos de teste descritos no capítulo de testes controlados.
4. Comparar os resultados com as classificações `APROVADO`, `FALHOU` e `NOTA`.

